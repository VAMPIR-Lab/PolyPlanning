Comment and uncomment the first function in nonsmooth_setup.jl to change sd function.  
Use this command to run experiment with 400 samples(4 maps * 100 initial conditions):  
include("./utilities/run_simple_experiment.jl")  

Result:  
Success rate for new sd function is 63.75% , 51% for old sd function.  
The warning "More derivatives than accounted for" is more likely to occur on new version.  
The case that mcp is solved but only one vertex of ego touches the obstacle is more likely to occur on new version.  
