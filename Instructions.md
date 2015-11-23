# Notes

3T data coregisteration: shift of the small images over the 60100 (slFOV)
maskshft is made from MRA data 3T mat file MaskMRACut  which has been shifted 1 to the right using shiftmask code.
maskshft is made from MRA data 7T mat file MaskMRACut  which has been shifted 0 to the left and 5 up using shiftmask code.
maybe MaskMRACut for 7T needs to be shrinked one pixel? used bwmorph( thin) in code MaskShrink and creating maskshftShrink7T
now ShiftMaskBack can be used to creat maskVX6071Shrink7T for example that is needed for unwrapping and Ensight generation.

## 3T
MRAMaskCut3T data-including slices: [19 60]

4040  data-including slices: [9 40(theortically)]


6040 [1+93 120+93]

4040 [1+93 120+93]

6071 [1+45 216+45]

## 7T
MRAMaskCut7T data-including slices: [16 60]

4040  data-including slices: [6 40(theortically)]


6040 [1+92 120+92]

4040 [1+92 120+92]

6071 [1+44 216+44]

ZeroPadding in both directions for both 3T and 7T: 
6040 padded 92
4040 padded 92 
6071 padded 44 


