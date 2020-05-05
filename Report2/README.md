This Readme file Shall guide you to Running the report codes.

# 1- Requirements:
Results can be found in:
* MultiVars_Comparison_results.txt for multivariable optimization comparison.
* SingleVars_Comparison_results.txt for single variable optimization comparison.

# 2- To reproduce the results, you may run:
* Compare_Multi_Vars.m : to reproduce  multivariable optimization comparison results. <br>
In this file you will find the implementation of the following algorithms: Fletcher-Reeves Conjugate Gradient Method, Marquardt Method, and Broyden-Fletcher-Goldfarb-Shanno algorithm.
* Compare_Single_Vars.m : to reproduce  singlevariable optimization comparison results. <br>
In this file you will find implementation of the following algorithms: Fibonnaci method, golden section method, quadratic interpolation method, cubic interpolation method.

# 3- Testing:
The two mentioned  codes are user-friendly. 
Search for "EDIT PART BELOW TO TEST" in any of the files mentioned in #2 and replace the objecative function there with your function to test on any other function.
Make sure the function is symbolic and that the symbols are defined (using syms var1 var2 ... varn)
