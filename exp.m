close all; clear all; clc;

syms x t g l m
eqn = (1-t*x)/(g-(1-t*x)*l) + t*x*(12*m-x*l)/(8*m*(m-x*l)) == 0
[solx, params, conds] = solve(eqn, x, 'ReturnConditions', true)

solx
params
conds
