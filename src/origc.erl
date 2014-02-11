%
% 2013-11-24 karlll <karl@ninjacontrol.com>
%

%
% Finding origc
%

-module(origc).

-compile([export_all]).


%% -------------------------------------------------------------------------- %%
%% Frequency                                                                  %%
%% -------------------------------------------------------------------------- %%

%
% Print most frequent k-mers of a specified length in a providede string
%
print_most_frequent_kmers([],_) ->
	io:format("None.");
print_most_frequent_kmers(_,0) ->
	io:format("None.");
print_most_frequent_kmers(String,Length) ->
	TopKmers = most_frequent_kmers(String,Length),
	{_,Freq,_} = hd(TopKmers),
	io:format("Most frequent K-mer(s), length ~p :~n",[Length]),
	lists:keymap(fun(El) -> io:format("~s~n",[atom_to_list(El)]) end, 1, TopKmers),
	io:format("Frequency : ~p~n",[Freq]).




%
% Returns a tuple list of the most frequent kmers of lengt Length in 
% the provided string, of the format
% [{kmer_atom, frequency}]
%
most_frequent_kmers(String,Length) ->
	Freqs = get_kmer_freq(String,Length),
	{_,TopFreq,_} = hd(Freqs),
	lists:filter(fun(El) -> case El of 
							{_,TopFreq,_} ->
								true;
							_ ->
								false
							end
				end,
				Freqs
				).

%
% Returns a list of the format
% [{kmer_atom,{frequency,position_list}}], sorted in descending frequency
%
get_kmer_freq(String,Length) ->
	get_kmer_freq(String,length(String),Length,0,[]).
get_kmer_freq([],_,_,_,Acc) ->
	lists:reverse(lists:keysort(2,Acc));
get_kmer_freq(String,StringLength,KmerLength,Pos,Acc) ->
	case StringLength of
		L when L < KmerLength ->
			lists:reverse(lists:keysort(2,Acc));
		_L ->
			{Kmer,_} = lists:split(KmerLength,String),
			Key = list_to_atom(Kmer),
			NewAcc = case lists:keyfind(Key,1,Acc) of
				false -> [{Key,1,[Pos]}|Acc]; % first entry
				{Key,Count,PosList} -> lists:keyreplace(Key,1,Acc,{Key,Count+1,[Pos|PosList]}) 
			end,
			get_kmer_freq(lists:nthtail(1,String),StringLength-1,KmerLength,Pos+1,NewAcc)
	end.


%% -------------------------------------------------------------------------- %%
%% Complement                                                                 %%
%% -------------------------------------------------------------------------- %%

%
% compl(S) -> S'
% Returns the complement of string S
%
compl($A) -> $T;
compl($T) -> $A;
compl($G) -> $C;
compl($C) -> $G;

compl(String) ->
	compl(String,[]).

compl([N|Tail],Acc) ->
	compl(Tail,[compl(N)|Acc]);
compl([],Acc) ->
	Acc.

%% -------------------------------------------------------------------------- %%
%% Find substrings                                                            %%
%% -------------------------------------------------------------------------- %%


print_substring_positions(SubString,String) ->
	Pos = find_substring(SubString,String),
	io:format("Positions :~n"),
	lists:foreach(fun(El) -> io:format("~p ",[El]) end, Pos).

%
% Find the positions of substring SubString in String
%
find_substring(SubString,String) ->
	LenSub = length(SubString),
	LenStr = length(String),
	lists:reverse(find_substring(SubString,String,LenSub,LenStr,0,[])).


find_substring(_,[],_,0,_,Acc) ->
	Acc; % end of string
find_substring(SubString,String,LenSub,LenStr,Pos,Acc) ->
	case LenStr of
		L when L < LenSub ->
			Acc; % no more possible matches
		_ ->	
			case lists:prefix(SubString,String) of
				true -> NewAcc = [Pos|Acc];
				false -> NewAcc = Acc
			end,
			find_substring(SubString,lists:nthtail(1,String),LenSub,LenStr-1,Pos+1,NewAcc)
	end.

%% -------------------------------------------------------------------------- %%
%% Find clumps                                                                %%
%% -------------------------------------------------------------------------- %%

%
% TODO: this is not very efficient, re-write!
%

print_find_clumps(String,WindowLength,KmerLength,ClumpNum) ->
	Clumps = find_clumps(String,WindowLength,KmerLength,ClumpNum),
	ClumpsStr = lists:map(fun(El) -> {K,_,_} = El, atom_to_list(K) end, lists:flatten(Clumps)),
	lists:foreach(fun(El) -> io:format("~s~n",[El]) end, lists:usort(ClumpsStr)).


%
% Find the k-mers of length KmerLength within a substring of WindowLength
% in String occuring ClumpNum number of times
%

find_clumps(String,WindowLength,KmerLength,ClumpNum) ->
	find_clumps(String,length(String),WindowLength,KmerLength,ClumpNum,0,[]).

find_clumps([],_,_,_,_,_,Acc) ->
	Acc;
find_clumps(String,StringLen,WindowLength,KmerLength,ClumpNum,OffsetCounter,Acc) ->
	case StringLen of
		L when L < WindowLength ->
			Acc;
		_L ->
			{Window,_} = lists:split(WindowLength,String),
			Clumps = clumps(Window,KmerLength,ClumpNum),
			NewAcc = case Clumps of
						[] -> Acc;
						_ -> [Clumps|Acc]
					 end,
			find_clumps(lists:nthtail(1,String),StringLen-1,WindowLength,KmerLength,ClumpNum,OffsetCounter+1,NewAcc)
	end.


clumps(String,KmerLength,ClumpNum) ->
	Freqs = get_kmer_freq(String,KmerLength),
	lists:filter(fun(El) -> case El of 
							{_,ClumpNum,_} -> true;
							_ -> false
						end
					end, Freqs).

%% -------------------------------------------------------------------------- %%
%% Skew                                                                       %%
%% -------------------------------------------------------------------------- %%

print_skew(String) ->
	S = skew(String),
	lists:foreach(fun(El) -> io:format("~p ", [El]) end, S).

%
% skew(S) -> Skew::list().
%
% Caluculate skew for string S
%

skew(String) when is_list(String) ->
	skew(String,[0]);
skew($A) -> 0;
skew($T) -> 0;
skew($C) -> -1;
skew($G) -> 1.

skew([],Acc) ->
	lists:reverse(Acc);
skew([H|Tail],Acc) ->
	skew(Tail,[skew(H)+hd(Acc)|Acc]).

print_min_skew(String) ->
	S = min_skew(String),
	lists:foreach(fun(El) -> io:format("~p ", [El]) end, S).

%
% min_skew(S) -> Pos::list().
%
% Find position in S which minimizes skew
%

min_skew(String) ->
	S = skew(String),
	Min = lists:min(S),
	min_skew(S,Min,0,[]).

min_skew([],_,_,Acc) ->
	lists:reverse(Acc);
min_skew([Min|Tail],Min,Count,Acc) ->
	min_skew(Tail,Min,Count+1,[Count|Acc]);
min_skew([_H|Tail],Min,Count,Acc) ->
	min_skew(Tail,Min,Count+1,Acc).

%% -------------------------------------------------------------------------- %%
%% Approx. pattern matching                                                   %%
%% -------------------------------------------------------------------------- %%

print_approx_match(Pattern,String,Mismatches) ->
	ApproxMatches = approx_match(Pattern,String,Mismatches),
	lists:foreach(fun(El) -> {Pos,_,_} = El, io:format("~p ", [Pos]) end, ApproxMatches).

write_approx_match_positions(Pattern,String,Mismatches) ->
	ApproxMatches = approx_match(Pattern,String,Mismatches),
	ResultStr = lists:map(fun(El) -> {Pos,_,_} = El, io_lib:fwrite("~p ", [Pos]) end, ApproxMatches),
	util:write_result(lists:flatten(ResultStr),"out.result").

%
% Find all patterns in string String which has at most Mismatches differences
% from Pattern.
%
% approx_match(P,S,M) -> [{Pos::integer(),Count::integer(),Pattern2::string()}]
%
approx_match(Pattern,String,Mismatches) ->
	LenPattern = length(Pattern),
	LenString = length(String),
	approx_match(Pattern,String,Mismatches,LenPattern,LenString,0,[]).

approx_match(_,[],_,_,_,_,Acc) ->
	lists:reverse(Acc);
approx_match(Pattern,String,Mismatches,LenPattern,LenString,Pos,Acc) ->
	case LenString of
		L when L < LenPattern ->
			lists:reverse(Acc);
		_ -> 
			{Pattern2,_} = lists:split(LenPattern,String),
			NewAcc = case count_mismatch(Pattern,Pattern2) of
						C when C =< Mismatches ->
							[{Pos,C,Pattern2}|Acc];
						_ ->
							Acc
					end,
			approx_match(Pattern,lists:nthtail(1,String),Mismatches,LenPattern,LenString-1,Pos+1,NewAcc)
	end.


count_mismatch(String1,String2) when length(String1) == length(String2) ->
	count_mismatch(String1,String2,0).

count_mismatch([],[],Count) ->
	Count;
count_mismatch([H|Tail1],[H|Tail2],Count) ->
	count_mismatch(Tail1,Tail2,Count);
count_mismatch([_H|Tail1],[_N|Tail2],Count)  ->
	count_mismatch(Tail1,Tail2,Count+1).









