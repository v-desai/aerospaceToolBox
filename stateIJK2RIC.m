function [stateRIC,covRIC] = stateIJK2RIC(refstate,state,cov)
% state [6x1] pos and vel
% cov [3x3] pos
rvec = state(1:3,1);
vvec = state(4:6,1);
Beta = IJK2RIC(rvec,vvec);
stateRIC = Beta*rvec;
covRIC = Beta*cov*Beta';