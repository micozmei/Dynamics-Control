function [ M ] = Rotation_Matrix( angle )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
M = [cos(angle), -sin(angle)
    sin(angle), cos(angle)];
end

