This solves, using Matlab, the ode system (6) presented in (v1 of) the following paper: https://arxiv.org/abs/2505.10622
See endnote [14]
Something like figure 3 will be made if this program is run as given. 
I package the initial guess functions chebfuns at a couple of points, so you'll need to install that matlab add on.
It is probably not needed, so you can rewrite that bitof the code to avoid, making them just ordinary functions that are input into the guess function/subroutine.
That is all.
