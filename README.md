Current status:

```
1. Being able to get a result for all the 9 cars which satisfied all the distance constraints between any two cars.

2. The computation time is very long. It takes approximately 30s to finish the iterations and get a converged results.
   Possible reasons:
   i: The trajectory variables' dimension is too high. In this case the dimension is 630.
   ii: There are too many constraints. For every time step, 9 cars will induce 9*8/2 = 36 distance constraints. And as there are
   in total 35 trajectory points, there are 35 * 36 = 1260 constraints, which severely slow down the optimization.

```

