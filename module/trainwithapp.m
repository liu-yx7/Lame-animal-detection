fm1=features;
[fm1, PS]=mapminmax(fm1');
fm1=fm1';
fm1=array2table(fm1);
fm3= cell2table(target);
fm2= cell2table(featureslabels);
fm=[fm1 fm2 fm3];
save ('fm');