
%
% 2013-12-03 karlll <karl@ninjacontrol.com>
%


-module(util).

-compile([export_all]).

%% -------------------------------------------------------------------------- %%
%% Utility functions                                                          %%
%% -------------------------------------------------------------------------- %%

write_result(Result,Filename) ->
	file:write_file(Filename, io_lib:fwrite("~p.\n", [Result])).

write_result_text(Result,Filename) ->
	file:write_file(Filename, io_lib:fwrite("~s\n", [Result])).

load_string(FileName) ->
	V = load_binary_string(FileName),
	L = binary_to_list(V),
	string:strip(L,both,$\n).

load_binary_string(FileName) ->
	{ok, BinaryString} = file:read_file(FileName),
	BinaryString.
	





get_input(Filename) ->
	{ok, F} = file:open(Filename,[read]),
	L = read_input(F),
	get_input(L,[]).

get_input(["Input:\n"|T],Acc) ->
	get_input(T,Acc);
get_input(["Input\n"|T],Acc) ->
	get_input(T,Acc);
get_input(["Output:\n"|_T],Acc) ->
	Acc;
get_input(["Output\n"|_T],Acc) ->
	Acc;
get_input([],Acc) ->
	Acc;
get_input([L|T], Acc) ->
	get_input(T,[string:strip(L,right,$\n)|Acc]).

read_input(File) ->
    case file:read_line(File) of
        eof        -> [];
        {ok, Line} -> [Line | read_input(File)]
    end.



%
% matrix, M x N (Row-Col)
% TODO: implement using stdlib array instead


new_mx(M,N) ->
	new_mx(M,N,0).
% new MxN matrix, init all values to InitValue
new_mx(M,N,InitValue) ->
	[ lists:duplicate(N,InitValue) || _X <- lists:seq(1,M) ].


% Set value a Row M, Col N in Rowlist to Value. Zero based.
mx(M,N,Value,Rowlist) ->
	mx(M,N,Value,0,Rowlist,[]).
mx(M,N,Value,RowC,[Row|T],Acc) when M =/= RowC ->
	mx(M,N,Value,RowC+1,T,[Row|Acc]);
mx(M,N,Value,RowC,[Row|T],Acc) when M =:= RowC ->
	NewRow = lists:sublist(Row,N) ++ [Value] ++ lists:nthtail(N+1,Row),
	lists:reverse(Acc) ++ [NewRow] ++ T.

%
% Set column N to column vector ColVector in matrix Mx. Zero based.
%
mx_col(N,ColVector,Mx) ->
	mx_col(N,0,ColVector,Mx).
mx_col(_,_,[],Acc) ->
	Acc;
mx_col(N,Count,[Value|ColT],Acc) ->
	mx_col(N,Count+1,ColT,mx(Count,N,Value,Acc)).


%
% Set row M to row vector RowVector in matrix Mx. Zero based. 
%
mx_row(M,RowVector,Mx) ->
	lists:sublist(Mx,M) ++ [RowVector] ++ lists:nthtail(M+1,Mx).

% return value in Nnth col, Mth row
% Row & col numbering is zero based. 
mx_at(M,N,RowList) ->
	%io:format("mx_at: M = ~p, N = ~p~n",[M,N]),
	lists:nth(N+1,lists:nth(M+1,RowList)).

% return Nth column vector in matrix. Zero based.
mx_col_at(N,RowList) ->
	lists:map(fun(Row) -> lists:nth(N+1,Row) end,RowList).

% return Mth row vector in matrix. Zero based.
mx_row_at(M,RowList) ->
	lists:nth(M+1,RowList).

% Return maximum element along with its position in matrix Mx; {Value,I,J}
mx_max(Mx) ->
	MxN = lists:zip(Mx,lists:seq(0,length(Mx)-1)),
	Max = lists:map(fun(M)->{R,I} = M,{MaxV,J} = mx_max_vector(R),{MaxV,I,J} end,MxN),
	Res = lists:last(lists:keysort(1,Max)),
	Res.

% Get max value and its position in vector; {MaxValue,Pos}
mx_max_vector(V) ->
	mx_max_vector(V,0,{-1,-1}).

mx_max_vector([],_I,Max) ->
	case Max of
		{-1,-1} ->
			none;
		_ ->
			Max
	end;

mx_max_vector([H|T],I,Max) ->
	{MaxV,_} = Max,
	NewMax = case H of
		H when H > MaxV ->
			{H,I};
		_ ->
			Max
	end,
	mx_max_vector(T,I+1,NewMax).


%
% max
%

v_max(A,B) when A > B -> A;
v_max(A,B) when B > A -> B;
v_max(A,B) when A =:= B -> A.

v_max(A,B,C) ->
	v_max(v_max(A,B),C).



print_path(Nodes) ->
	string:join(Nodes,"->").

% find position of first occurence of El in List.
idx(El,List) ->
	idx(El,List,1).
idx(El,[El|_T],Idx) ->
	Idx;
idx(El,[_H|T],Idx) ->
	idx(El,T,Idx+1);
idx(_El,[],_Idx) ->
	not_found.

% pop list until head is element El
pop_until(El,[El|T]) ->
	[El|T];
pop_until(_El,[]) ->
	[];
pop_until(El,[H|T]) ->
	pop_until(El,T).


% Extract a sublist starting from first occurence of element StartEl
% to first occurrence of element EndEl. Precond: StartEl must preceed EndEl
sublist2(StartEl,EndEl,List) ->
	StartIdx = idx(StartEl,List),
	EndIdx = idx(EndEl,List),
	lists:sublist(List,StartIdx,EndIdx-StartIdx+1).












