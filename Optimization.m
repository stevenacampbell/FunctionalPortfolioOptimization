%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Optimization of Portfolio Generating Function %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% IMPORT DATA %%%

%Import (u,v) data where v=u_kr_k/<u,r>
%The data is of dimention nxT where n is the number of stocks and T
%is the number of periods
filename0='u_data.csv';
vecu = readmatrix(filename0);
filename1='v_data.csv';
vecv = readmatrix(filename1);

%Import the grid points for the piecewise affine function
filename2='grid.csv';
g = readmatrix(filename2);

%Get grid increment
dg = diff(g);

%%% PROBLEM SETTINGS %%%

% Note: The objective function optimized here is rescaled by a factor of T
% and does not include the term corresponding to a re-weighted diversity.
% See the associated paper for the original and most general statement of the problem.
    
%Option for the Monotonicity of Port. Weights 1-Yes, 0-No
Mon = 0; 

%Set the value at which the function must be 0
%This specifies a constraint and the value should lie in the grid
val=0.5;

%Set beta for optimization
beta = 1000000;

%Recover problem parameters from imported data
dims = size(vecu); %Get input dimensions
T = dims(2); %# of Periods
n = dims(1); %# of Stocks
d=length(g); %Dimension of Grid

%%% OPTIMIZATION %%%

%Specify solver
cvx_solver Mosek_3;

%Formulate optimization problem (specify objective and constraints)
cvx_begin
    variable x(d)
    expression s(T)
    expression dx(d-1)
    dx = x(2:d)-x(1:d-1);
    %Build objective function
    for i=1:T
        s(i)=1;
        dvu=(vecv(:,i)-vecu(:,i));
        for j=1:n
            k=MatchGrid(vecu(j,i),g);
            if k==d
               k = d-1;
            end
            s(i)=s(i)+(dvu(j)/(dg(k)*n))*dx(k);
        end
    end
    obj = sum(log(s)); 
    %Specify problem type (i.e. Maximize objective function)
    maximize(obj)
    %Build constraints
    subject to
        for i=1:(d-2)
            %Exponential Concavity Constraints
            w=dg(i)/(dg(i+1)+dg(i));
            -dx(i)+log(w*exp(dx(i+1)+dx(i))+(1-w))<=0;
            %Beta-Smoothness Constraints
            abs((dg(i)/dg(i+1))*dx(i+1)-dx(i))<=beta*dg(i)^2;  
        end
        %Endpoint Constraints
        (dx(1))^2-beta*dg(1)^2<=0; 
        (dx(d-1))^2-beta*dg(d-1)^2<=0;
        %Fixed Function Value Constraint
        idx=MatchGrid(val,g);
        x(idx)==0; 
        %Monotonicity Constraint
        if(Mon==1) 
            dx(2)>=0;
            for i=2:(d-2)
                (g(i+1)/g(i))*(dg(i)/dg(i+1))*dx(i+1)-dx(i)>=0;
            end
        end
cvx_end

%%% PLOT AND OUTPUT SOLUTION %%%

plot(g,x)

filename = 'solution.csv';
writematrix(x, filename) 


%%% HELPER FUNCTION %%%

%Function to match point in [0,1] to a gridpoint
%Here the convention is to choose the closest gridpoint but the point can
%also be matched according to the closest gridpoint less than it, etc.
function idx = MatchGrid(p,grid)
    idx = find(grid>=p,1);
    if idx>1 && abs(grid(idx-1)-p)<=abs(grid(idx)-p)
        idx = idx-1;
    end
end

