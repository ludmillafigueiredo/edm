module pollinators

#using

#export

function createPollinator(x,y)
	global pollinator_id_counter =  Int64(1)

	Pollinator(
	id = string(pollinator_id_counter),
	sp = 'sp1',
	genotype = 'AA',
	#TODO traits depend on genotype
	biomass = 0.01,
	disp_par = (1,1),
	#TODO should the pollinator object hold an interaction matrix or verify it when polinizing? (the latter, but how)
	poll_access = 'flower',
	location = (x,y)
	)
end

function pollinator_growth(pollinator::Pollinator)
	grown_mass += growthrate * pollinator.biomass^(3/4) * exp(-Ea/(Boltz*fragment[pollinator.location].temperature)) #TODO same growthrate for insects and plants?
	#TODO can this reference to the fragment be a problem? how will it know which fragment (only  a problem if using the 3D array arrangement for landscape, actually)

	pollinator.biomass += grown_mass
end

function pollinator_dispersal() #TODO make an interaction function and use differnt methods depending on partner?
	# pollinator carries copy of gamet of visited flower

end


end
