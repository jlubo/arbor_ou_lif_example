#!/bin/python3

import arbor
import numpy as np
import time
import matplotlib.pyplot as plt

#####################################
# TestRecipe
# Implementation of Arbor simulation recipe
class TestRecipe(arbor.recipe):

	# constructor
	def __init__(self):
		arbor.recipe.__init__(self)
		
		self.props = arbor.neuron_cable_properties() # initialize the cell properties
		cat = arbor.load_catalogue("./custom-catalogue.so") # load the catalogue of custom mechanisms
		cat.extend(arbor.default_catalogue(), "") # add the default catalogue
		self.props.catalogue = cat

		self.runtime = 20000 # runtime of the simulation in ms of biological time
		self.dt = 0.2 # duration of one timestep in ms

		self.stim1_prot = {"time_start": 10,
			               "scheme": "TRIPLET"} # protocol for stimulation, "time_start": starting time, "scheme": protocol type
		self.stim2_prot = {"time_start": 15,
			               "scheme": "ONEPULSE"} # protocol for stimulation, "time_start": starting time, "scheme": protocol type
		self.bg_prot = {"time_start": 0,
			            "scheme": "FULL"} # protocol for background input, "time_start": starting time, "scheme": protocol type

	# cell_kind
	# Defines the kind of the neuron given by gid
	# - gid: global identifier of the cell
	# - return: type of the cell
	def cell_kind(self, gid):
		
		return arbor.cell_kind.cable # note: implementation of arbor.cell_kind.lif is not ready to use yet

	# cell_description
	# Defines the morphology, cell mechanism, etc. of the neuron given by gid
	# - gid: global identifier of the cell
	# - return: description of the cell
	def cell_description(self, gid):

		assert gid == 0 # only one neuron in use

		# cylinder morphology
		tree = arbor.segment_tree()
		radius = 1e-10 # radius of cylinder (in µm)
		height = 2*radius # height of cylinder (in µm)
		tree.append(arbor.mnpos,
		            arbor.mpoint(-height/2, 0, 0, radius),
		            arbor.mpoint(height/2, 0, 0, radius),
		            tag=1)
		labels = arbor.label_dict({"center": "(location 0 0.5)"})
		area_m2 = 2 * np.pi * (radius * 1e-6) * (height * 1e-6) # surface area of the cylinder in m^2 (excluding the circle-shaped ends, since Arbor does not consider current flux there)
		area_cm2 = 2 * np.pi * (radius * 1e-4) * (height * 1e-4) # surface area of the cylinder in cm^2 (excluding the circle-shaped ends, since Arbor does not consider current flux there)
		i_factor = (1e-9/1e-3) / area_cm2 # conversion factor from nA to mA/cm^2; for point neurons
		C_mem = 1e-9 # neuronal capacitance in F
		c_mem = C_mem / area_m2 # specific capacitance in F/m^2, computed from absolute capacitance of point neuron
		V_rev = -65.0 # reversal potential in mV
		V_th = -55.0 # spiking threshold in mV
		
		# cell mechanism
		decor = arbor.decor()
		decor.set_property(Vm=V_rev, cm=c_mem)
		mech_neuron = arbor.mechanism("lif")
		h_0 = 4.20075 # initial synaptic weight in nC
		R_leak = 10.0 # leak resistance in MOhm
		tau_mem = R_leak*10**9 * C_mem # membrane time constant in ms
		tau_syn = 5.0 # synaptic time constant in ms
		mech_neuron.set("R_leak", R_leak)
		mech_neuron.set("R_reset", 1e-10)
		mech_neuron.set("I_0", 0) # set to zero (background input is applied via OU process ou_bg)
		mech_neuron.set("i_factor", i_factor)
		mech_neuron.set("V_rev", V_rev)
		mech_neuron.set("V_reset", -70.0)
		mech_neuron.set("V_th", V_th)
		mech_neuron.set("t_ref", 2.0) # refractory time in ms
		decor.paint('(all)', arbor.density(mech_neuron))

		print("area =", area_m2, "m^2")
		print("i_factor =", i_factor, "(mA/cm^2) / (nA)")
		print("c_mem =", c_mem, "F/m^2")
		print("tau_mem =", tau_mem, "ms")

		# background input current to all neurons (described by Ornstein-Uhlenbeck process)
		mech_ou_bg = arbor.mechanism('ornstein_uhlenbeck_gen')
		mech_ou_bg.set('mean', 0.15) # mean current in nA
		mech_ou_bg.set('stdev', 0.05) # standard deviation in nA*s^(1/2)
		mech_ou_bg.set('tau', tau_syn) # synaptic time constant in ms
		decor.place('"center"', arbor.synapse(mech_ou_bg), "ou_bg")

		# input current described by Ornstein-Uhlenbeck process accounting for a population of neurons
		mech_ou_stim = arbor.mechanism('ornstein_uhlenbeck_pop')
		mech_ou_stim.set('N', 25) # number of neurons in the putative population
		mech_ou_stim.set('f', 100) # firing rate of the putative neurons in Hz 
		mech_ou_stim.set('w_out', h_0/R_leak) # synaptic weight in nC
		mech_ou_stim.set('tau', tau_syn) # synaptic time constant in ms
		decor.place('"center"', arbor.synapse(mech_ou_stim), "ou_stim")		

		# place spike detector
		decor.place('"center"', arbor.spike_detector(V_th), "spike_detector")
			
		return arbor.cable_cell(tree, labels, decor)
		
	# connections_on
	# Defines the list of incoming synaptic connections to the neuron given by gid
	# - gid: global identifier of the cell
	# - return: connections to the given neuron
	def connections_on(self, gid):
		
		return []
	
	# getOUStimulusGenerators
	# Returns arbor.event_generator instances that describe the specifics of a given
	# input/stimulation protocol for a mechanism implementing an Ornstein-Uhlenbeck process. 
	# Here, if the value of the 'weight' parameter is 1, stimulation is switched on, 
	# whereas if it is -1, stimulation is switched off.
	# - protocol: protocol that is being used for stimulation
	# - label: label of the mechanism that receives the the stimulation
	# - return: a list of arbor.event_generator instance of the stimulation protocol
	def getOUStimulusGenerators(self, protocol, label):

		prot_name = protocol['scheme'] # name of the protocol (defining its structure)
		start_time = protocol['time_start'] # time at which the stimulus starts in s

		if prot_name == "ONEPULSE":
			# create regular schedules to implement a stimulation pulse that lasts for 0.1 s
			stim_on = arbor.event_generator(
			          label,
			          1,
			          arbor.regular_schedule(start_time*1000, self.dt, start_time*1000 + self.dt))
			stim_off = arbor.event_generator(
			           label,
			           -1,
			           arbor.regular_schedule((start_time + 0.1)*1000, self.dt, (start_time + 0.1)*1000 + self.dt))  
			
			return [stim_on, stim_off]
			
		elif prot_name == "TRIPLET":
			# create regular schedules to implement pulses that last for 0.1 s each
			stim1_on = arbor.event_generator(
			           label,
			           1,
			           arbor.regular_schedule(start_time*1000, self.dt, start_time*1000 + self.dt))
			stim1_off = arbor.event_generator(
			            label,
			            -1,
			            arbor.regular_schedule((start_time + 0.1)*1000, self.dt, (start_time + 0.1)*1000 + self.dt))
			stim2_on = arbor.event_generator(
			           label,
			           1,
			           arbor.regular_schedule((start_time + 0.5)*1000, self.dt, (start_time + 0.5)*1000 + self.dt))
			stim2_off = arbor.event_generator(
			            label,
			            -1,
			            arbor.regular_schedule((start_time + 0.6)*1000, self.dt, (start_time + 0.6)*1000 + self.dt))
			stim3_on = arbor.event_generator(
			           label,
			           1,
			           arbor.regular_schedule((start_time + 1.0)*1000, self.dt, (start_time + 1.0)*1000 + self.dt))
			stim3_off = arbor.event_generator(
			            label,
			            -1,
			            arbor.regular_schedule((start_time + 1.1)*1000, self.dt, (start_time + 1.1)*1000 + self.dt))
				  
			return [stim1_on, stim1_off, stim2_on, stim2_off, stim3_on, stim3_off]
			
		elif prot_name == "FULL":
			# create a regular schedule that lasts for the full runtime
			stim_on = arbor.event_generator(
			           label,
			           1,
			           arbor.regular_schedule(start_time*1000, self.dt, start_time*1000 + self.dt))
				  
			return [stim_on]

		return []

	# event_generators
	# Event generators for input to synapses
	# - gid: global identifier of the cell
	# - return: events generated from Arbor schedule
	def event_generators(self, gid):
		inputs = []

		# OU stimulus input #1
		stim = self.getOUStimulusGenerators(self.stim1_prot, "ou_stim")
		inputs.extend(stim)		

		# OU stimulus input #2
		stim = self.getOUStimulusGenerators(self.stim2_prot, "ou_stim")
		inputs.extend(stim)

		# OU background input
		stim = self.getOUStimulusGenerators(self.bg_prot, "ou_bg")
		inputs.extend(stim)
			
		return inputs
		
	# global_properties
	# Sets properties that will be applied to all neurons of the specified kind
	# - gid: global identifier of the cell
	# - return: the cell properties 
	def global_properties(self, kind): 

		assert kind == arbor.cell_kind.cable # assert that all neurons are technically cable cells

		return self.props
	
	# num_cells
	# - return: the total number of cells in the network
	def num_cells(self):
		
		return 1

	# probes
	# - gid: global identifier of the cell
	# - return: the probes on the given cell
	def probes(self, gid):

		# probe membrane potential, total current, and external input currents
		return [arbor.cable_probe_membrane_voltage('"center"'), \
				arbor.cable_probe_total_ion_current_cell(), \
				arbor.cable_probe_point_state_cell("ornstein_uhlenbeck_pop", "I_ou"), \
				arbor.cable_probe_point_state_cell("ornstein_uhlenbeck_gen", "I_ou")]

#####################################
if __name__ == '__main__':

	#####################################
	# set up and run simulation
	recipe = TestRecipe()

	random_seed = int(time.time()*10000) # get random seed from system clock
	print("random_seed = " + str(random_seed))
	
	alloc = arbor.proc_allocation(threads=1, gpu_id=None) # select one thread and no GPU (default; cf. https://docs.arbor-sim.org/en/v0.7/python/hardware.html#arbor.proc_allocation)
	context = arbor.context(alloc, mpi=None) # constructs a local context without MPI connection
	domains = arbor.partition_load_balance(recipe, context) # constructs a domain_decomposition that distributes the cells in the model described by an arbor.recipe over the distributed and local hardware resources described by an arbor.context (cf. https://docs.arbor-sim.org/en/v0.5.2/python/domdec.html#arbor.partition_load_balance)
	sim = arbor.simulation(recipe, context, domains, seed = random_seed)

	reg_sched = arbor.regular_schedule(0, recipe.dt, recipe.runtime) # create schedule for recording 

	# set handles to probe membrane potential and currents
	gid = 0
	handle_mem = sim.sample((gid, 0), reg_sched) # membrane potential
	handle_tot_curr = sim.sample((gid, 1), reg_sched) # total current
	handle_stim_curr = sim.sample((gid, 2), reg_sched) # stimulus current
	handle_bg_curr = sim.sample((gid, 3), reg_sched) # background current

	sim.record(arbor.spike_recording.all)
	sim.run(tfinal=recipe.runtime, dt=recipe.dt)

	#####################################
	# get traces and spikes from simulator
	loc = 0 # for single-compartment neurons, there is only one location

	data_mem, data_curr = [], []

	data, _ = sim.samples(handle_mem)[loc]
	times = data[:, 0]

	# read out neuron data
	data_mem.append(data[:, 1])

	# get total membrane current
	data_tot_curr, _ = sim.samples(handle_tot_curr)[loc]
	data_curr.append(np.negative(data_tot_curr[:, 1]))

	# get Ornstein-Uhlenbeck currents
	for hnd in [handle_stim_curr, handle_bg_curr]:
		A = None # array to append to data_curr
		if len(sim.samples(hnd)) > 0:
			data_ou_curr, _ = sim.samples(hnd)[loc]
			if len(data_ou_curr[0]) > 0:
				A = data_ou_curr[:, 1]
		if A is None:
			if hnd == handle_stim_curr:
				print("No samples for stimulus current.")
			elif hnd == handle_bg_curr:
				print("No samples for background current.")
			A = np.zeros(len(data[:, 0]))
		data_curr.append(A)

	# get spikes
	spike_times = sim.spikes()['time']

	#####################################
	# assemble, store, and plot data
	data_header = "Time, V, I_stim, I_bg"
	data_stacked = np.column_stack([times, \
		                            np.transpose(data_mem), np.transpose(data_curr)])

	np.savetxt("traces.txt", data_stacked, fmt="%.4f", header=data_header)
	np.savetxt("spikes.txt", spike_times, fmt="%.4f")

	fig, axes = plt.subplots(nrows=4, ncols=1, sharex=False, figsize=(10, 10)) # create figure with 'num_rows' subfigures
	axes[0].set_title("random_seed = " + str(random_seed))
	axes[0].plot(data_stacked[:,0], data_stacked[:,1], label='V', color='C3')
	axes[0].plot(spike_times, -55*np.ones(len(spike_times)), '.', color='blue', markersize=1)
	axes[1].plot(data_stacked[:,0], data_stacked[:,2], label='I_tot', color='C1')
	axes[2].plot(data_stacked[:,0], data_stacked[:,3], label='I_stim', color='C2')
	axes[3].plot(data_stacked[:,0], data_stacked[:,4], label='I_bg', color='C0')
	axes[3].set_xlabel("Time (ms)")
	axes[0].set_ylabel("V (mV)")
	axes[1].set_ylabel("I_tot (nA)")
	axes[2].set_ylabel("I_stim (nA)")
	axes[3].set_ylabel("I_bg (nA)")
	fig.savefig("traces.svg")
	#fig.savefig("traces.png", dpi=800)
	#plt.show()
	
