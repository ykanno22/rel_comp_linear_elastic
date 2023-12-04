Copyright (C) Yoshihiro Kanno, 2021. All rights reserved.
===============================================================================
Y. Kanno: "Computation-with-confidence for static elasticity: data-driven approach with order statistics." ZAMM, Vol.103, e202100482 (2023).

The codes can run on Matlab ver.9.7 with CVX ver. 2.2 and SeDuMi ver. 1.3.4. 
Also, the following Matlab codes available at 
   <http://au.mathworks.com/support/books/book118765.html>
are required for FEM. 
  - elem_q4_fun.m
  - fmlin.m
  - formbee.m
  - formdsig.m
  - gauss.m
  - Q4_mesh_fun.m
These FEM codes are contained in the following book:
  - A. Khennane:
    Introduction to Finite Element Analysis Using MATLAB and Abaqus.
    CRC Press, Boca Raton (2013).

To reproduce the result presented in section 5.2.2, run the following sequentially: 
  - lin_orthotropic_data_set.m
  - lin_6d_eig_loop.m
Part of these codes and the functions used by them was based on [Khennane, 2013] cited above.
