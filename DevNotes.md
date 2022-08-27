# Notes about development

These notes are a couple of attention point for developper. Reviewing them can help you if you are looking for a bug. 

## R8Ki or DbKi?
It is often unclear to me where to use `R8Ki`or `DbKi` for the code. Bad usage can result in some seg fault appearing randomly due to memory compromission. 
In particular, in order to cast `C_DOUBLE` into `R8Ki` safely, it is adviced to do change the type using the function `real(..,Type)` AND to use the array numbers : `x(1:3)`

```f90
real(c_double),intent(in)  :: input_c(3)
real(R8KI)                 :: output_f(3)
output_f(1:3) = real(NacOri_C(1:3),R8Ki)
```