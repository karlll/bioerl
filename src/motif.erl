%
% 2013-11-27 karlll <karl@ninjacontrol.com>
%

-module(motif).
-compile([export_all]).




%% -------------------------------------------------------------------------- %%
%% Motif enum                                                                 %%
%% -------------------------------------------------------------------------- %%

%
% NB: I don't think motif_enum works as it is supposed to...
%
print_motif_enum(Dna,K,D) ->
    Me = motif_enum(Dna,K,D),
    lists:foreach(fun(M) -> io:format("~s ",[M]) end, Me).

motif_enum(Dna,K,D) ->
	Kmers = get_kmers_from_list(Dna,K),
    motif_enum(Kmers,Dna,D,[]).

motif_enum([],_Dna,_D,Acc) ->
    lists:sort(sets:to_list(sets:from_list(lists:concat(Acc)))); % flatten and remove duplicates
motif_enum([Kmer|Ktail],Dna,D,Acc) ->
    MutKmers = get_mutations(Kmer,D), % Kmers differing from Kmer by at most D mutations
    NewAcc = find_candidates(MutKmers,Dna,D),
    motif_enum(Ktail,Dna,D,[NewAcc|Acc]).

find_candidates(Kmers,Dna,D) ->
    find_candidates(Kmers,Dna,D,[]).

find_candidates([],_Dna,_D,Acc) ->
    Acc;

find_candidates([Kmer|Ktail],Dna,D,Acc) ->
    MutKmers = get_mutations(Kmer,D),
    NewAcc = case is_present(MutKmers,Dna) of
            true ->
                [Kmer | Acc];
            _ ->
                Acc
            end,
    find_candidates(Ktail,Dna,D,NewAcc).

is_present(Kmers,Dna) ->
    lists:all(
            fun(DnaL) -> 
                
                lists:any(
                    fun(Kmer) ->
                        case string:str(DnaL,Kmer) of % is Kmer a substring of Dnal
                            0 -> false;
                            _ -> true
                        end
                    end,
                    Kmers)
            end,
            Dna).

get_kmers_from_list(LString,K) ->
    get_kmers_from_list(LString,K,[]).

get_kmers_from_list([],_K,Acc) ->
    Acc; 
get_kmers_from_list([L|LTail],K,Acc) ->
    get_kmers_from_list(LTail,K,get_kmers(L,K)++Acc).


% get unique kmers of length K in string String
get_kmers(String,K) ->
	get_kmers(String,K,[]).

get_kmers(String,K,Acc) when length(String) < K ->
	sets:to_list(sets:from_list(Acc)); % unique kmers
get_kmers(String,K,Acc) when length(String) >= K ->
	{Kmer,_} = lists:split(K,String),
	get_kmers(lists:nthtail(1,String),K,[Kmer|Acc]).

% get all combination of Kmer with at most N (0..N) mutated Bps
get_mutations(Kmer,N) ->
	get_mutations(Kmer,N,[]).

get_mutations(_Kmer,0,Acc) ->
    sets:to_list(sets:from_list(Acc)); % remove duplicates
get_mutations(Kmer,N,Acc) ->
    get_mutations(Kmer,N-1,mutate(Kmer,N,"ATGC") ++ Acc).    

% Get all strings based on String with N positions mutated with values from Values list

mutate(String,N,Values) ->
    P = comb(N,lists:seq(1,length(String))),
    B = comb_rep(N,Values),
    {_,MTuples,_} = sofs:product(sofs:set(P),sofs:set(B)),
    lists:map(fun(T) -> 
            {Pos,Vals} = T,
            replaceat(String,Pos,Vals)
        end,
        MTuples).

% replace elerments in string at positions Positions with corresponding
% values from Values
replaceat(String,[P|Ptail], [V|Vtail]) ->
    replaceat(replaceat(String,P,V),Ptail,Vtail);

replaceat(String,[],[]) ->
    String;

replaceat(String,Position,Value) when is_integer(Position), is_integer(Value) ->
    replaceat(String,1,Position,Value,[]).


replaceat([H|Tail],Pos,Pos,Value,Acc) ->
    lists:reverse(Acc) ++ [Value|Tail];

replaceat([H|Tail],Pos,NPos,Value,Acc) ->
    replaceat(Tail,Pos+1,NPos,Value,[H|Acc]).    

% Get all combinations of size N from List ( snippet from Rosetta Code http://rosettacode.org/wiki/Combinations#Erlang)
comb(0,_) ->
    [[]];
comb(_,[]) ->
    [];
comb(N,[H|T]) ->
    [[H|L] || L <- comb(N-1,T)]++comb(N,T).

comb_rep(0,_) ->
    [[]];
comb_rep(_,[]) ->
    [];
comb_rep(N,[H|T]=S) ->
    [[H|L] || L <- comb_rep(N-1,S)]++comb_rep(N,T).

perms([]) -> [[]];
perms(L)  -> [[H|T] || H <- L, T <- perms(L--[H])].


kmer_perms_rep(K) ->
    lists:map(fun(El) -> string:right(dec_to_bp(El),K,$A) end, % pad with 'A'
        lists:seq(0,trunc(math:pow(4,K))-1)).



dec_to_bp(Int) ->
    dec_to_bp(integer_to_list(Int,4),[]).

dec_to_bp([],Acc) ->
    lists:reverse(Acc);
dec_to_bp([H|T],Acc) ->
    dec_to_bp(T,[bp(H)|Acc]).

bp($0) -> $A;
bp($1) -> $T;
bp($2) -> $G;
bp($3) -> $C.

%% -------------------------------------------------------------------------- %%
%% Mean string                                                                %%
%% -------------------------------------------------------------------------- %%

% find mininum Hamming distance in string Dna for pattern Pattern
d(Pattern, Dna) ->
    Kmers = get_kmers(Dna,length(Pattern)),
    lists:min(lists:map(fun(Kmer) -> hamming(Pattern,Kmer) end, Kmers)).



% find mininum Hamming distance in string Dna for pattern Pattern, retrun tuple
% {Distance,Kmer}
d2(Pattern, Dna) ->
    Kmers = get_kmers(Dna,length(Pattern)),
    AllDK = lists:keysort(1,lists:map(fun(Kmer) -> {hamming(Pattern,Kmer),Kmer} end, Kmers)),
    {MinD,_} = hd(AllDK),
    lists:filter(fun(El) -> 
                    {D,_} = El, 
                    case D of
                        MinD ->
                            true;
                        _ -> 
                            false
                    end
                end,
                AllDK).

d2sum(Pattern,DnaL) ->
    lists:foldl(fun(El,Acc) -> {D,_} = hd(d2(Pattern,El)), Acc+D end, 0, DnaL).


% calculate hamming distance between P1 & P2
hamming(P1,P2) when length(P1) == length(P2) ->
    hamming(P1,P2,0).

hamming([],[],Acc) ->
    Acc;
hamming([P|T1],[P|T2],Acc) ->
    hamming(T1,T2,Acc);
hamming([_P|T1],[_Q|T2],Acc) ->
    hamming(T1,T2,Acc+1).


print_median_string(Dna,K) ->
    MedianStr = median_string(Dna,K),
    lists:foreach(fun(El) -> io:format("~s~n",[El]) end, MedianStr).

median_string(Dna,K) ->
    Kmers = kmer_perms_rep(K),
    median_string(lists:nthtail(1,Kmers),Dna,d2sum(hd(Kmers),Dna),[hd(Kmers)]).

median_string([],_Dna,_BestD,Acc) ->
    Acc;

median_string([Kmer|Ktail],Dna,BestD,Acc) ->
    case d2sum(Kmer,Dna) of
        X when X < BestD ->
            median_string(Ktail,Dna,X,[Kmer]);
        Y when Y > BestD ->
            median_string(Ktail,Dna,BestD,Acc);
        Z when Z =:= BestD ->
            median_string(Ktail,Dna,Z,[Kmer|Acc])
    end.


%% -------------------------------------------------------------------------- %%
%% Profile most probable k-mer                                                %%
%% -------------------------------------------------------------------------- %%

print_prob_kmer(String,K,PMatrix,ColMap) ->
    HighestProb = prof_prob_kmer(String,K,PMatrix,ColMap),
    lists:foreach(fun(El) -> io:format("~s~n",[El]) end, HighestProb).

prof_prob_kmer(String,K,PMatrix,ColMap) when K =:= length(PMatrix) ->
    Kmers = get_kmers(String,K),
    Probs = lists:map(fun(Kmer) -> {get_kmer_prob(Kmer,PMatrix,ColMap),Kmer} end,Kmers),
    SortedProbs = lists:reverse(lists:keysort(1,Probs)),
    {HighestP,_} = hd(SortedProbs),
    lists:filtermap(fun(P) -> 
                        {N,Km} = P, 
                        case N of
                            N when N =:= HighestP ->
                                {true,Km};
                            _ ->
                                false
                        end
                    end, SortedProbs).



get_kmer_prob(Kmer,PMatrix,ColMap) ->
    %io:format("Kmer = ~s ",[Kmer]),
    get_kmer_prob(Kmer,1,PMatrix,ColMap,1).
get_kmer_prob([],_N,_PMatrix,_ColMap,Acc) ->
    %io:format(" = ~p~n",[Acc]),
    Acc;
get_kmer_prob([K|Ktail],N,PMatrix,ColMap,Acc) ->
    P = pmatrix(K,N,ColMap,PMatrix),
    %io:format("~p ",[P]),
    get_kmer_prob(Ktail,N+1,PMatrix,ColMap,Acc * P).
            


% get probability for base B appearing at position N in a kmer of length length(Rows)
% ColMap maps bases to cols [{$A,1},{$T,2},{$G,3},{$C,4}]
pmatrix(B,N,ColMap,Rows) ->
    Row = lists:nth(N,Rows),
    ColN = proplists:get_value(B,ColMap,none),
    lists:nth(ColN,Row).


def_colmap() -> [{$A,1},{$C,2},{$G,3},{$T,4}].


