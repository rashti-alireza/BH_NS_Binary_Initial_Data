#########################################
## some tips for setting the param files:
#########################################

o. Set params "BH_radius_update_weight" and "BH_spin_update_weight" 
   to zero for resolution 20x20x20 to stop changing the inner boundary 
   condition on the apparent horizon. Otherwise, it might cause instability. 
   However, if you need more precision on BH mass and spin, you may set 
   these weight params to a smaller value like 0.1.

o. For high spinning BH cases increase the iterations, my suggestion is:
   12(x200)->14(x100)->16(x100)->18(x100)->20(x50)->22(x50)

o. Memory requirements on SplitCubedSpherical(BH+NS) grid where BH is
   excised:
   . resolution 22x22x22 ~ 100 GB.
   . resolution 24x24x24 ~ 120 GB.

