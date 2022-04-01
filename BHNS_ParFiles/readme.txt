#########################################
## some tips for setting the param files:
#########################################


## conformally flat BHNS:

o. A safe choice of iteration for a production run can be:

BH_radius_update_weight = 0.(x10)->0.3(x390)->0.
BH_spin_update_weight   = 0.(x10)->0.3(x390)->0.
n_a = 12(x200)->14(x100)->16(x100)->18(x100)->20(x70)->22(x50)
n_b = 12(x200)->14(x100)->16(x100)->18(x100)->20(x70)->22(x50)
n_c = 12(x200)->14(x100)->16(x100)->18(x100)->20(x70)->22(x50)

However, this number of iteration is an over-killing for 
simple cases of BHNS, i.e., zero spin, well resolved BH and 
not too close objects. These cases can be treated with 
a fewer number of iteration; see the other param files 
for some examples.

o.  Memory requirements on SplitCubedSpherical(BH+NS) grid where BH is
excised and UMFPACK is used:
. resolution 22x22x22 ~ 100 GB.
. resolution 24x24x24 ~ 120 GB.

o.  Set params "BH_radius_update_weight" and "BH_spin_update_weight" 
to zero for resolution 18x18x18 to stop changing the inner boundary 
condition on the apparent horizon. Otherwise, it might cause instability. 
However, if you need more precision on BH mass and spin, you may set 
these weight params to a smaller value like 0.1.

o.  To drive the P_ADM to zero for difficult cases, such as high spinning BH
and NS in arbitrary directions, it helps to increase the number of
iteration at high resolutions. I have also tried to increase the param
"BHNS_P_ADM_control_update_weight" but it did not help and from previous
experiments it might cause instability. Thus, keep the weight the same but
increase the iterations.

o. Note: smaller size BHs, for example, high spin BHs or small mass ones,
required more resolutions to drive the constraints small enough.


## conformally curved BHNS:

