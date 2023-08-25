module runtests

include(test_dir*"matrix_IO_test.jl")
include(test_dir*"moments_test.jl")

include(test_dir*"cp_matrices_test.jl")
include(test_dir*"cp_model_test.jl")

include(test_dir*"nn_matrices_test.jl")
include(test_dir*"nn_model_test.jl")

end




