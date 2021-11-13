classdef DirectSolver < Solver

    methods (Static, Access = public)
        function solucio = solve(RHS, LHS)
            solucio = RHS \ LHS;
        end
    end

end