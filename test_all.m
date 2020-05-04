% test script for all tests

p=Person(0.0);
p.ConsistencyCheck;
DNA=3*rand(1,p.GetNumberOfActions*p.GetNumberOfStates);
DNA=round(DNA);   

disp('test_person')
tic;test_person;toc
tic;test_person(DNA);toc

disp('test_community')
tic;test_community;toc
tic;test_community(100);toc
tic;test_community(100,20);toc
tic;test_community(100,20,DNA);toc


disp('test_stoc_community')
tic;test_stoc_community;toc
tic;test_stoc_community(DNA);toc