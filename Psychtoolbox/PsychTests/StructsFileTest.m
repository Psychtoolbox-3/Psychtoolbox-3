% StructsFileTest
%
% Test the routines for reading and writing numeric structs
% into and out of text files.
%
% 6/16/03  dhb  Wrote it.

% Clear
clear all; close all;

% Fill in some structs
theStructs(1).field1 = 'joe';
theStructs(1).field2 = 2;
theStructs(1).field3 = 3;
theStructs(1).field4 = 4;

theStructs(2).field1 = 'jim';
theStructs(2).field2 = 6;
theStructs(2).field3 = 7;
theStructs(2).field4 = 8;

theStructs(3).field1 = 'bob';
theStructs(3).field2 = 10;
theStructs(3).field3 = 11;
theStructs(3).field4 = 12;

% Write the array out
WriteStructsToText('foo.txt',theStructs);

% Read it back in
newStructs = ReadStructsFromText('foo.txt');

% Delete test file.
delete foo.txt

