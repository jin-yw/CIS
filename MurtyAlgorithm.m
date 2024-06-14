function [ costs, xSet ] = MurtyAlgorithm( C, k )
% The implementation of Murty's Algorithm with Hungairan Algorithm
% Reference:
% M. L. Miller, H. S. Stone and I. J. Cox, "Optimizing Murty's ranked assignment method," in IEEE Transactions on Aerospace and Electronic Systems, vol. 33, no. 3, pp. 851-862, July 1997, doi: 10.1109/7.599256.
% 1. Find the best solution, S0, to P0.
[cost0,x0] = HungarianAlgorithm(C);
if k == 1
    % If only 1 best assignment is required:
    costs = cost0;
    xSet = cell(1,1);
    xSet{1} = x0;
else
    % 2. Initialize a priority queue of problem/solution pairs to contain
    % only <P0,S0>. The top pair on this queue will always be the pair with
    % the lowest-cost solution.
    Data.C = C;
    Data.S = x0;
    Data.cost = cost0;
    q0 = dlnode(Data);
    start = q0;
    % 3. Clear the list of solutions to be retured.
    costs = [];
    xSet = cell(0,0);
    % 4. For i = 1 to k, or until the priority queue of problem/solution
    % pairs is empty:
    i = 1;
    while i <= k && ~isempty(q0)
        % 4.1 Take the top problem/solution pair, <P,S>, off the queue.
        TopPair = start.Data;
        start = q0.Next;
        removeNode(q0);
        q0 = start;
        % 4.2 Add S to the list of solutions to be returned.
        if ~isempty(costs)
            if TopPair.cost == costs(end)
                continue;
            else
                costs(i) = TopPair.cost;
                xSet{i} = TopPair.S;
            end
        else
            costs(1) = TopPair.cost;
            xSet{1} = TopPair.S;
        end
        % 4.3 For each triple, <y,z,l>, found in S:
        [idx1, idx2] = find(TopPair.S == 1);
        idxLen = length(idx1);
        for j = 1:idxLen
            % 4.3.1 Let P' = P
            tempC = TopPair.C;
            % 4.3.2 Remove the triple <y,z,l> from P'
            tempC(idx1(j), idx2(j)) = inf;
            % 4.3.3 Find the best solution, S', to P'.
            existFlag = true;
            [tempCost,tempX0] = HungarianAlgorithm(tempC);
            if isempty(tempX0)
                existFlag = false;
            end
%             try
%                 [tempCost,tempX0] = HungarianAlgorithm(tempC);
%             catch
%                 existFlag = false;
%             end
            % 4.3.4 If S' exists:
            if existFlag
                % 4.3.4.1 Add <P',S'> onto the queue of
                % problem/solution pairs.
                cruiser = start;
                exitFlag = false;
                tempData.C = tempC;
                tempData.S = tempX0;
                tempData.cost = tempCost;
                q1 = dlnode(tempData);
                while ~exitFlag
                    if isempty(cruiser)
                        cruiser = q1;
                        start = q1;
                        q0 = start;
                        exitFlag = true;
                    else
                        if tempData.cost < cruiser.Data.cost
                            insertBefore(q1,cruiser);
                            if isempty(q1.Prev)
                                start = q1;
                                q0 = q1;
                            end
                            exitFlag = true;
                        end
                        if ~isempty(cruiser.Next) && ~exitFlag
                            cruiser = cruiser.Next;
                        elseif isempty(cruiser.Next) && ~exitFlag
                            insertAfter(q1,cruiser);
                            exitFlag = true;
                        end
                    end
                end
            end
            % 4.3.5 From P, remove triples that include y, and all
            % triples that include z, except for <y,z,l> itself.
            tempValue = TopPair.C(idx1(j),idx2(j));
            TopPair.C(idx1(j),:) = inf;
            TopPair.C(:,idx2(j)) = inf;
            TopPair.C(idx1(j),idx2(j)) = tempValue;
        end
        i = i + 1;
    end
end

end

