
% This script performs fundamental social network analysis tasks on the
% collaboration graph for the ICMB conference for the time period between
% 2002 and 2013.
% Clear workspace and command window.
clc
clear all

% Set the period of years.
Years = [2002:1:2013];
YearsNum = length(Years);

% Load weight matrices for each year.
for year = Years
    filename = strcat(['ICMB-' num2str(year) '.mat']);
    load(filename);
end;

% Load authors' names.
load('authors.mat');

% Set a container for the weight matrices of all years.
ICMB = cell(1,numel(Years));

% ICMB{1} = array_2002;
% ICMB{2} = array_2003;
% ICMB{3} = array_2004;
% ICMB{4} = array_2005;
% ICMB{5} = array_2006;
% ICMB{6} = array_2007;
% ICMB{7} = array_2008;
% ICMB{8} = array_2009;
% ICMB{9} = array_2010;
% ICMB{10} = array_2011;
% ICMB{11} = array_2012;
% ICMB{12} = array_2013;

for y = 1:YearsNum

    ICMB{y} = eval(genvarname(strcat(['array_' num2str(Years(y))])));

end;

% Get the number of nodes N.
N = size(ICMB{1},1);

% Construct the overall graph weight matrix.
W = zeros(N,N);
for y = 1:1:YearsNum
    W = W + ICMB{y};
end;

% Set up a vector of indices pointing to the diagonal elements of the
% weight matrix W.
Idiag = [1:N+1:N*N];

% Re-initialize the overall weight matrix W so that fundamental social
% network analysis tasks can be performed. W should be a binary adjacency 
% matrix so that W[i,j] = 1 indicates the presence of an edge between 
% authors i and j. Moreover, the diagonal elements of W should also be set 
% to zero.
Wo = W;
Wo(Wo>1) = 1;
Wo(Idiag) = 0;

% Extract connected components of co-authorship network.
% Initially extract the corresponding NeighboursList.
NL = NeighboursList(Wo);
C = ConnectedComponents(NL);
N=10; % here we set the number of the top N authors to be returned
for i=1:length(C)
    X = ['For Component ',num2str(i),' :'];
    disp(X);
    if (length(C{i})<N)
        disp('Not enough nodes in component');
        continue;
    end;
    %intialize Wo matrix to the correct format and as input in FloydWarshall algorithm 
    W1=Wo(C{i},C{i}); 
    W1(W1==0) = Inf;
    W1(eye(size(W1))~=0) = 0;
    [D,P]=FloydWarshall(W1);
    D(~isfinite(D)) = 0; % set infinite distances to 0 to make more efficient calculations of centrality measures
    %display(num2str(size(W1,2)));
    closeness_centrality=sum(W1~=0,2)./sum(D,2); %estimate closeness centrality by dividing connected nodes with the total distance from all nodes
    eccentricity_centrality=1./max(D,[],2); % calculating eccentricity as 1 over the maximum distance of the shortest paths distance
    MeasureName = 'Closeness Centrality';
    MeasureValues = closeness_centrality;
    ReportTopNAuthors(MeasureValues,MeasureName,N,authors(C{i},:)); % print top N authors of component
    MeasureName = 'Eccentricity Centrality';
    MeasureValues = eccentricity_centrality;
    ReportTopNAuthors(MeasureValues,MeasureName,N,authors(C{i},:)); % print top N authors of component
    
    %betweenness centrality
    %Betweenness centrality quantifies the number of times a node
    %acts as a bridge along the shortest path between two other nodes
    %In order to find the betweenness centrality of node i
    %we need to calculate gjk which is the number of shortest paths between all node
    %pairs (j,k) and gjk(i) which is the number of shortest paths between node pairs (j,k)
    %going through i
    %then betweenness centrality equals to (gjk(i)/gjk)
    %number of nodes in component
    W2=Wo(C{i},C{i});
    n=length(W2);      
    %initialize the centrality
    I=eye(n)~=0;       
    %initialize shortest path length
    d=1; 
    %number of paths of length |d|
    NPd=W2;     
    %number of shortest paths of length |d|
    NSPd=NPd;         
    %number of shortest paths of any length
    NSP=NSPd; NSP(I)=1;  
    %length of shortest paths
    L=NSPd; L(I)=1;           	

    %calculate NSP and L
    %while number of shortest paths of length |d| contains 1
    while find(NSPd,1)
        % increase the path length by 1
        d=d+1;
        % calculate the product of the adjacency matrix multipied by the 
        % number of paths of length |d| matrix 
        NPd=NPd*W2;
        % element wise multiplication 
        NSPd=NPd.*(L==0);
        % increase number of shortest paths of any length by umber of
        % shortest paths of length |d| 
        NSP=NSP+NSPd;
        % increase the length of shortest path by path length
        L=L+d.*(NSPd~=0);
    end
    %distance length for disconnected vertices is inf
    L(~L)=inf; L(I)=0;         
    %number of shortest paths for disconnected vertices is 1
    NSP(~NSP)=1;                

    Gt=W2.';
    %initialize dependency matrix
    DP=zeros(n); 
    %graph diameter
    diam=d-1;                  	

    %calculate dependencies
    for d=diam:-1:2
        DPd1=(((L==d).*(1+DP)./NSP)*Gt).*((L==(d-1)).*NSP);
        %DPd1: dependencies on vertices |d-1| from source
        DP=DP + DPd1;           
    end
    betweenness_centrality = sum(DP,1); 
    %normalize: divide bc with the the highest possible value of
    %betweenness
    %The highest possible value of betweenness is ((n-2)*(n-1))/2
    %this is the max number node pairs in
    %an undirected graph except node if the node was directed highest
    %possible value would be (n-2)*(n-1)
    betweenness_centrality = 2*betweenness_centrality./((n-2)*(n-1));
    MeasureName = 'Betweenness Centrality';
    MeasureValues = betweenness_centrality;
    ReportTopNAuthors(MeasureValues,MeasureName,N,authors(C{i},:));  
    
    %eigenvector centrality
    %It is a measure of the influence of a node in a network.
    %it is calculated as The centrality of each vertex is proportional to the
    %sum of the centralities of its neighbors
    %Higher 1st eigenvector scores indicate more central nodes to the 
    %main pattern of distances among all of the nodes, lower values 
    %indicate that actors are more peripheral. 
    
    % number of nodes in this component
    n = length(W2);
    
    if n < 1000
        % eig() is faster but only appropriate for small size matrices which
        % fit well in memory. Eig does not work with sparse matrixes
        [V,D] = eig(W2);
    else
        % sparse()converts a full matrix into sparse form by 
        % squeezing out any zero elements. If a matrix contains many zeros, 
        % converting the matrix to sparse storage saves memory
        % returns a vector of W2 eigenvalues
        % D is a diagonal matrix of w2's eigenvalues
        % V is a full matrix whose columns are the corresponding
        % eigenvectors
        [V,D] = eigs(sparse(W2));
    end
    % ~ represents an output that is discarded
    % here we want to know the index of the maximum eigenvalue
    % not the actual value
    [~,idx] = max(diag(D));
    % keep the column of V which contains the eigenvectors coresponding to
    % the maximum eigenvalue
    ec = abs(V(:,idx));
    % reshape the eigenvector's matrix into a 1 column matrix
    % each row will now contain the eigenvector_centrality of the 
    % author with id equal to the row number of the eigenvector_centrality matrix 
    eigenvector_centrality = reshape(ec, length(ec), 1);
    MeasureName = 'Eigenvector Centrality';
    MeasureValues = eigenvector_centrality;
    ReportTopNAuthors(MeasureValues,MeasureName,N,authors(C{i},:));    
    
end;   
