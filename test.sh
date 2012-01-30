#! /bin/sh -e
#
echo "Running make test in cppad_codegen/test"
cd cppad_codegen/test 
make test
cd /home/joao/Development/CppADCodeGen
#
echo "Running make test in speed/cppad"
cd speed/cppad 
make test
cd /home/joao/Development/CppADCodeGen
#
echo "Running make test in speed/double"
cd speed/double 
make test
cd /home/joao/Development/CppADCodeGen
#
echo "Running make test in speed/example"
cd speed/example 
make test
cd /home/joao/Development/CppADCodeGen
#
echo "Running make test in speed/profile"
cd speed/profile 
make test
cd /home/joao/Development/CppADCodeGen
#
echo "Running make test in example"
cd example 
make test
cd /home/joao/Development/CppADCodeGen
#
echo "Running make test in introduction/get_started"
cd introduction/get_started 
make test
cd /home/joao/Development/CppADCodeGen
#
echo "Running make test in introduction/exp_apx"
cd introduction/exp_apx 
make test
cd /home/joao/Development/CppADCodeGen
#
echo "Running make test in multi_thread"
cd multi_thread 
make test
cd /home/joao/Development/CppADCodeGen
#
echo "Running make test in print_for"
cd print_for 
make test
cd /home/joao/Development/CppADCodeGen
#
echo "Running make test in test_more"
cd test_more 
make test
cd /home/joao/Development/CppADCodeGen
exit 0
