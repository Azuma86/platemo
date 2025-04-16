classdef originalRWMOP < PROBLEM
% <multi> <real> <constrained>
% A simple n-dimensional constrained problem with a parameter 'Alpha'.
%
% The dimension D can be arbitrary (default: 2). The domain is [0,1]^D.
% We introduce two constraints, both controlled by Alpha:
%  1) sum(x_i) <= Alpha
%  2) sum((x_i - 0.5)^2) <= (Alpha/2)^2
%
% You can adjust 'Alpha' to make the feasible region larger (bigger Alpha)
% or smaller (smaller Alpha).

    properties
        Alpha = 0.03;  % パラメータ (可行領域を制御)
    end

    methods
        %% Initialization
        function Setting(obj)
            % 目的は2つ (M=2)
            obj.M = 2;
            % Dが未指定の場合は2次元をデフォルトとする
            if isempty(obj.D)
                obj.D = 2;
            end
            % 設計変数の範囲は [0,1]^D
            obj.lower    = zeros(1,obj.D);
            obj.upper    = ones(1,obj.D);
            obj.encoding = ones(1,obj.D);
            obj.number = 0;
        end
        
        %% Evaluate multiple solutions
        function Population = Evaluation(obj,varargin)
            % 取り出し
            x = varargin{1};       % 各行が1個体 (次元D)
            N = size(x,1);         % 個体数

            % ------ 目的関数例 ------
            % f1: sum of squares
            % f2: sum of squares from 1
            % (今回「考えなくてよい」とあるので、簡易的に定義)
            f1 = sum(x.^2, 2);
            f2 = sum((x - 1).^2, 2);
            f  = [f1, f2];

            % ------ 制約の定義 ------
            alpha = obj.Alpha;
            g = zeros(N,2);

            % (1) sum of x_i <= alpha
            g(:,1) = sum(x,2) - alpha;

            % (2) sum_i ((x_i - 0.5)^2) <= (alpha/2)^2
            g(:,2) = sum((x - 0.5).^2, 2) - (alpha/2)^2;

            % 不等式制約: g(:,i) <= 0
            Population = SOLUTION(x, f, g, varargin{2:end});
            obj.FE = obj.FE + N;  % 評価回数を更新(PlatEMOの管理用)
        end
        
        %% Generate a point for hypervolume calculation (任意)
        function R = GetOptimum(obj,~)
            % ここでは単純に [0,0] を返しておく
            R = [0,0];
        end
    end
end