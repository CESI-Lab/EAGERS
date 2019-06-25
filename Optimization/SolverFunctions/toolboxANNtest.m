function [ANN, accTrain, accTest] = toolboxANNtest(depth, features_train, features_test, sols_train, sols_test)
%this function outputs a trained artificial neural network made from
%matlab's neural network toolbox
%inputs: depth: integer value for number of layers in network,
%features_train: matrix of normalized feature examples for training,
%features_test: matrix of normalized feature examples ofor testing,
%sols_train: vector of zeros and ones for two classes for training, 
%sols_test: vector of zeros and ones for two classes for testing
%outputs: ANN: network structure, accTrain: training accuracy, accTest:
%testing accuracy

%create the network
nF = length(features_train(1,:));%number of features
ANN = patternnet(depth);
[ANN, tr] = train(ANN, features_train, sols_train);
% %ANN = network(nF, depth);
% %[ANN, ~] = train(ANN, features_train, sols_train);
%classify using the trained network
trainSols = ANN(features_train);
%split into two classes
trainSols = (trainSols>0.5);
%calculate training accuracy 
accTrain = nnz(trainSols == sols_train)/length(sols_train);
%classify test set
testSols = ANN(features_test);
%split into two classes
testSols = (testSols>0.5);
%calculate test accuracy
accTest = nnz(testSols == sols_test)/length(sols_test);

end