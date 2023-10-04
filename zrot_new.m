function [M] = zrot(angle)

c = cos(pi*angle/180);
s = sin(pi*angle/180);

M = [c s 0; -s c 0; 0 0 1];