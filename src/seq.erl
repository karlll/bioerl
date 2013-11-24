-module(seq).

-compile([export_all]).

%% -------------------------------------------------------------------------- %%
%% Translation                                                                %%
%% -------------------------------------------------------------------------- %%

%
% Translate RNA-string StringR to amino acid string
%
% translate(RNA::string()) -> Peptide::string()  
%
translate(StringR) ->
	translate(StringR,[]).

translate([],Acc) ->
	lists:reverse(Acc);
translate([X,Y,Z|Tail],Acc) ->
	case codon(X,Y,Z) of
		stop ->
			lists:reverse(Acc);
		A ->
			translate(Tail,[A|Acc])
	end.

codon($A,$A,$A) -> $K;
codon($A,$A,$C) -> $N;
codon($A,$A,$G) -> $K;
codon($A,$A,$U) -> $N;
codon($A,$C,$A) -> $T;
codon($A,$C,$C) -> $T;
codon($A,$C,$G) -> $T;
codon($A,$C,$U) -> $T;
codon($A,$G,$A) -> $R;
codon($A,$G,$C) -> $S;
codon($A,$G,$G) -> $R;
codon($A,$G,$U) -> $S;
codon($A,$U,$A) -> $I;
codon($A,$U,$C) -> $I;
codon($A,$U,$G) -> $M;
codon($A,$U,$U) -> $I;
codon($C,$A,$A) -> $Q;
codon($C,$A,$C) -> $H;
codon($C,$A,$G) -> $Q;
codon($C,$A,$U) -> $H;
codon($C,$C,$A) -> $P;
codon($C,$C,$C) -> $P;
codon($C,$C,$G) -> $P;
codon($C,$C,$U) -> $P;
codon($C,$G,$A) -> $R;
codon($C,$G,$C) -> $R;
codon($C,$G,$G) -> $R;
codon($C,$G,$U) -> $R;
codon($C,$U,$A) -> $L;
codon($C,$U,$C) -> $L;
codon($C,$U,$G) -> $L;
codon($C,$U,$U) -> $L;
codon($G,$A,$A) -> $E;
codon($G,$A,$C) -> $D;
codon($G,$A,$G) -> $E;
codon($G,$A,$U) -> $D;
codon($G,$C,$A) -> $A;
codon($G,$C,$C) -> $A;
codon($G,$C,$G) -> $A;
codon($G,$C,$U) -> $A;
codon($G,$G,$A) -> $G;
codon($G,$G,$C) -> $G;
codon($G,$G,$G) -> $G;
codon($G,$G,$U) -> $G;
codon($G,$U,$A) -> $V;
codon($G,$U,$C) -> $V;
codon($G,$U,$G) -> $V;
codon($G,$U,$U) -> $V;
codon($U,$A,$A) -> stop;
codon($U,$A,$C) -> $Y;
codon($U,$A,$G) -> stop;
codon($U,$A,$U) -> $Y;
codon($U,$C,$A) -> $S;
codon($U,$C,$C) -> $S;
codon($U,$C,$G) -> $S;
codon($U,$C,$U) -> $S;
codon($U,$G,$A) -> stop;
codon($U,$G,$C) -> $C;
codon($U,$G,$G) -> $W;
codon($U,$G,$U) -> $C;
codon($U,$U,$A) -> $L;
codon($U,$U,$C) -> $F;
codon($U,$U,$G) -> $L;
codon($U,$U,$U) -> $F.






% How many possible RNA strings encode Peptide PString ?
count_possible_strings(PString) ->
	count_possible_strings(PString,1).

count_possible_strings([A|PTail],Count) ->
	count_possible_strings(PTail,Count*length(aa(A)));
count_possible_strings([],Count) ->
	Count.


% Amino acid to codons
% aa(A::char) -> [Codons::list()]

aa($K) -> [[$A,$A,$A],[$A,$A,$G]];
aa($N) -> [[$A,$A,$C],[$A,$A,$U]];
aa($T) -> [[$A,$C,$A],[$A,$C,$C],[$A,$C,$G],[$A,$C,$U]];
aa($R) -> [[$A,$G,$A],[$A,$G,$G],[$C,$G,$A],[$C,$G,$C],[$C,$G,$G],[$C,$G,$U]];
aa($S) -> [[$A,$G,$C],[$A,$G,$U],[$U,$C,$A],[$U,$C,$C],[$U,$C,$G],[$U,$C,$U]];
aa($I) -> [[$A,$U,$A],[$A,$U,$C],[$A,$U,$U]];
aa($M) -> [[$A,$U,$G]];
aa($Q) -> [[$C,$A,$A],[$C,$A,$G]];
aa($H) -> [[$C,$A,$C],[$C,$A,$U]];
aa($P) -> [[$C,$C,$A],[$C,$C,$C],[$C,$C,$G],[$C,$C,$U] ];
aa($L) -> [[$C,$U,$A],[$C,$U,$C],[$C,$U,$G],[$C,$U,$U],[$U,$U,$A],[$U,$U,$G]];
aa($E) -> [[$G,$A,$A],[$G,$A,$G]];
aa($D) -> [[$G,$A,$C],[$G,$A,$U]];
aa($A) -> [[$G,$C,$A],[$G,$C,$C],[$G,$C,$G],[$G,$C,$U]];
aa($G) -> [[$G,$G,$A],[$G,$G,$C],[$G,$G,$G],[$G,$G,$U]];
aa($V) -> [[$G,$U,$A],[$G,$U,$C],[$G,$U,$G],[$G,$U,$U]];
aa($Y) -> [[$U,$A,$C],[$U,$A,$U]];
aa($C) -> [[$U,$G,$C],[$U,$G,$U]];
aa($W) -> [[$U,$G,$G]];
aa($F) -> [[$U,$U,$C],[$U,$U,$U]];
aa(stop) -> [[$U,$A,$A],[$U,$A,$G],[$U,$G,$A]].




%% -------------------------------------------------------------------------- %%
%% Aa naming functions                                                        %%
%% -------------------------------------------------------------------------- %%


aa_3_to_1(Name3) ->
	aa_3_to_1(Name3,[]).

aa_3_to_1([A,B,C|[]],Acc) when is_integer(C) ->
	lists:reverse([aa_name_to_1([A,B,C])|Acc]);
aa_3_to_1([A,B,C,$-|Tail],Acc) ->
	aa_3_to_1(Tail,[aa_name_to_1([A,B,C])|Acc]);
aa_3_to_1([],Acc) ->
	lists:reverse(Acc).

% Aa 3 letter abbrev. to single letter code
aa_name_to_1(Name3) when is_list(Name3), length(Name3) =:= 3 ->
	case lists:keyfind(Name3,2,aa_name()) of
		{_,Name3,Name1} ->
			Name1;
		false ->
			invalid_name
	end;

% Aa full name to single letter code
aa_name_to_1(Name) when is_list(Name) ->
	case lists:keyfind(Name,1,aa_name()) of
		{Name,_,Name1} ->
			Name1;
		false ->
			invalid_name
	end.

% Aa full name to 3 letter abbrev.
aa_name_to_3(Name) when is_list(Name) ->
	case lists:keyfind(Name,1,aa_name()) of
		{Name,Name3,_} ->
			Name3;
		false ->
			invalid_name
	end;

% Aa single letter code to 3 letter abbrev.
aa_name_to_3(Name1) when is_integer(Name1) ->
	case lists:keyfind(Name1,3,aa_name()) of
		{_,Name3,Name1} ->
			Name3;
		false ->
			invalid_name
	end.

% Aa 3 letter abbrev to full name.
aa_full_name(Name3) when is_list(Name3) and length(Name3) == 3 ->
	case lists:keyfind(Name3,2,aa_name()) of
		{Name,Name3,_} ->
			Name;
		false ->
			invalid_name
	end;

% Aa single letter code to full name
aa_full_name(Name1) when is_integer(Name1) ->
	case lists:keyfind(Name1,3,aa_name()) of
		{Name,_,Name1} ->
			Name;
		false ->
			invalid_name
	end.



aa_name() -> 
	[{"Alanine","Ala",$A},
	{"Arginine","Arg",$R},
	{"Asparagine","Asn",$N},
	{"Aspartic acid","Asp",$D},
	{"Cysteine","Cys",$C},
	{"Glutamic acid","Glu",$E},
	{"Glutamine","Gln",$Q},
	{"Glycine","Gly",$G},
	{"Histidine","His",$H},
	{"Isoleucine","Ile",$I},
	{"Leucine","Leu",$L},
	{"Lysine","Lys",$K},
	{"Methionine","Met",$M},
	{"Phenylalanine","Phe",$F},
	{"Proline","Pro",$P},
	{"Serine","Ser",$S},
	{"Threonine","Thr",$T},
	{"Tryptophan","Trp",$W},
	{"Tyrosine","Tyr",$Y},
	{"Valine","Val",$V},
	{"Selenocysteine","Sec",$U},
	{"Pyrrolysine","Pyl",$O},
	{"Asparagine or aspartic acid","Asx",$B},
	{"Glutamine or glutamic acid","Glx",$Z},
	{"Leucine or Isoleucine","Xle",$J},
	{"Unspecified or unknown amino acid","Xaa",$X}].

%% -------------------------------------------------------------------------- %%
%% Peptide naming functions                                                   %%
%% -------------------------------------------------------------------------- %%

peptide_mass_list(Peptides) ->
	peptide_mass_list(Peptides,[]).

peptide_mass_list([],Acc) ->
	Acc;

peptide_mass_list([P|Tail],Acc) ->
	peptide_mass_list(Tail,[peptide_mass_string(P)|Acc]).

peptide_mass_list_uniq(Peptides) ->
	Str = peptide_mass_list(Peptides),
	sets:to_list(sets:from_list(Str)).



%
% Convert a peptide of amino acid characters to a string of integers corresponding
% to the mass of the amino acid, separated by "-"
%
peptide_mass_string(PString) ->
	peptide_mass_string(PString,[]).

peptide_mass_string([],Acc) ->
	lists:flatten(Acc);

% this is ug-ley:
peptide_mass_string([H|[]],Acc) when is_integer(H) ->
	peptide_mass_string([],Acc ++ [integer_to_list(mass(H))]);
peptide_mass_string([H,N|[]],Acc) when is_integer(H), is_integer(N) ->
	peptide_mass_string([],Acc ++ [integer_to_list(mass(H)),$-,integer_to_list(mass(N))]);
peptide_mass_string([H,N|Tail],Acc) when is_integer(H), is_integer(N) ->
	peptide_mass_string(Tail,Acc ++ [integer_to_list(mass(H)),$-,integer_to_list(mass(N)),$-]).



%% -------------------------------------------------------------------------- %%
%% Peptide encoding                                                           %%
%% -------------------------------------------------------------------------- %%

rna(String) ->
	rna(String,[]).

rna([],Acc) ->
	lists:reverse(Acc);
rna([$T|Tail],Acc) ->
	rna(Tail,[$U|Acc]);
rna([H|Tail],Acc) ->
	rna(Tail,[H|Acc]).


print_p_enc(String,PString) ->
	Enc = p_enc(String,PString),
	lists:foreach(fun(El) -> io:format("~s~n",[El]) end, lists:reverse(Enc)).

%
% Return the substrings (and reverse complements) in String that encodes peptide PString
%

p_enc(String, PString) ->
	LS = length(String),
	LP = length(PString),
	p_enc(String,PString,LS,LP,[]).

p_enc([],_,_,_,Acc) ->
	lists:reverse(Acc);
p_enc(String,PString,LenString,LenPString,Acc) ->
	case LenString of
		L when L < LenPString * 3 -> % codon len = 3
			Acc;
		_ ->
			{Substr,_} = lists:split(LenPString * 3,String),
			SubstrCompl = origc:compl(Substr),
			DoesEnc = fun(S1,S2,T) -> 
						{r_enc_p(S1,T),r_enc_p(S2,T)} 
					  end,
			NewAcc = case DoesEnc(rna(Substr),rna(SubstrCompl),PString) of 
						{true,false} ->
							[Substr|Acc];
						{false,true} ->
							[Substr|Acc];
						_ ->
							Acc
					end,
			p_enc(lists:nthtail(1,String),PString,LenString -1, LenPString, NewAcc)
	end.		 	



r_enc_p([],[]) ->
	true;
r_enc_p([A,B,C|Rtail],[R|Ptail]) when is_integer(A),is_integer(B), is_integer(C), is_integer(R) ->
	case codon(A,B,C) of
		R ->
			r_enc_p(Rtail,Ptail);
		_ ->
			false
	end;

r_enc_p(_X,_Y) ->
	false.
	
%% -------------------------------------------------------------------------- %%
%%  Mass calc                                                				  %%
%% -------------------------------------------------------------------------- %%


mass(String) when is_list(String) ->
	mass(String,0);

mass($G) -> 57;
mass($A) -> 71;
mass($S) -> 87;
mass($P) -> 97;
mass($V) -> 99;
mass($T) -> 101;
mass($C) -> 103;
mass($I) -> 113;
mass($L) -> 113;
mass($N) -> 114;
mass($D) -> 115;
mass($K) -> 128;
mass($Q) -> 128;
mass($E) -> 129;
mass($M) -> 131;
mass($H) -> 137;
mass($F) -> 147;
mass($R) -> 156;
mass($Y) -> 163;
mass($W) -> 186.

mass([],Acc) ->
	Acc;
mass([H|Tail],Acc) when is_integer(H)->
	mass(Tail,Acc+mass(H)).


print_cyclo_spectrum(PString) ->
	CS = cyclo_spectrum(PString),
	lists:foreach(fun(El) -> io:format("~p ",[El]) end, CS).
%
% Calculate the theoretical spectrum of a cyclical peptide String
% [ 0 , m1, m2, .. mn, mass(P) ]
% InclPeptide == true includes peptide in list
%
cyclo_spectrum(PString) ->
	cyclo_spectrum(PString,false).

cyclo_spectrum(PString,InclPeptide) ->
	case InclPeptide of
		true -> cyclo_spectrum(cyclic_perm(PString),[{PString,mass(PString)}],true);
		_ -> cyclo_spectrum(cyclic_perm(PString),[mass(PString)],InclPeptide)
	end.

cyclo_spectrum([H|Tail],Acc,true) ->
	cyclo_spectrum(Tail,[{H, mass(H)}|Acc],true);

cyclo_spectrum([H|Tail],Acc,InclPeptide) ->
	cyclo_spectrum(Tail,[mass(H)|Acc],InclPeptide);

cyclo_spectrum([],Acc,true) ->
	lists:keysort(2,[{zero,0} | Acc]);
cyclo_spectrum([],Acc,_InclPeptide) ->
	lists:sort([0 | Acc]).


%
% Calculate the theoretical spectrum of a linear peptide String
% InclPeptide == true includes peptide in list
%

lin_spectrum(PString) ->
	lin_spectrum(PString,false).

lin_spectrum(PString,InclPeptide) ->
	case InclPeptide of
		true -> lin_spectrum(linear_perm(PString),[{PString,mass(PString)}],true);
		_ -> lin_spectrum(linear_perm(PString),[mass(PString)],InclPeptide)
	end.

lin_spectrum([H|Tail],Acc,true) ->
	lin_spectrum(Tail,[{H, mass(H)}|Acc],true);

lin_spectrum([H|Tail],Acc,InclPeptide) ->
	lin_spectrum(Tail,[mass(H)|Acc],InclPeptide);

lin_spectrum([],Acc,true) ->
	lists:keysort(2,[{zero,0} | Acc]);
lin_spectrum([],Acc,_InclPeptide) ->
	lists:sort([0 | Acc]).


%
% Get all ordered substrings, cyclic search
%

cyclic_perm(String) ->
	cyclic_perm(String,length(String),[]).

cyclic_perm(_,0,Acc) ->
	Acc;
	
cyclic_perm(String,Pos,Acc) ->
	Str = get_subst_wrap(String,Pos,length(String)-1),
	cyclic_perm(String,Pos-1,Str ++ Acc).


linear_perm(String) ->
	linear_perm(String,length(String),[]).

linear_perm(_,0,Acc) ->
	Acc;

linear_perm(String,Pos,Acc) ->
	Str = get_subst(String,Pos),
	linear_perm(String,Pos-1,Str ++ Acc).

%
% Get all substring from String starting at pos StartPos
%
get_subst(String,StartPos) ->
	get_subst(String,StartPos,1,length(String),[]).

get_subst(String,Pos,Len,MaxLen,Acc) when Len+Pos > MaxLen+1 ->
	lists:reverse(lists:delete(String,Acc)); % we will not include the whole string
get_subst(String,Pos,Len,MaxLen,Acc) when Pos+Len =< MaxLen+1 ->
	get_subst(String,Pos,Len+1,MaxLen,[lists:sublist(String,Pos,Len) | Acc]).

%
% Get all substrings up to length MaxLen from position StartPos.
% The search is done from left to right, and wraps from end of string to the beginning.
% 

get_subst_wrap(String,StartPos,MaxLen) ->
	get_subst_wrap(String,StartPos,1,MaxLen,[]).

get_subst_wrap(_String,_Pos,Len,MaxLen,Acc) when Len > MaxLen ->
	Acc;

get_subst_wrap(String,Pos,Len,MaxLen,Acc) when Len =< MaxLen, Pos+Len =< length(String)  ->  % no wrap
	get_subst_wrap(String,Pos,Len+1,MaxLen,[ lists:sublist(String,Pos,Len) | Acc ]);

get_subst_wrap(String,Pos,Len,MaxLen,Acc) when Len =< MaxLen, Pos+Len > length(String)  ->  % wrap
	LenA = length(String) - (Pos-1),
	LenB = Len - LenA,
	A = lists:sublist(String,Pos,LenA),
	B = lists:sublist(String,1,LenB),
	get_subst_wrap(String,Pos,Len+1,MaxLen,[ A ++ B | Acc ]).







%
% Cyclopeptide sequencing
%

print_cyclopeptide_seq(Spec) ->
	Candidates = cyclopeptide_seq(Spec),
	lists:foreach(fun(El) -> io:format("~s ",[El]) end, peptide_mass_list_uniq(Candidates)).


cyclopeptide_seq(Spec) ->
	cyclopeptide_seq(Spec,expand_peptide(""),[]).


cyclopeptide_seq(Spec,List,Acc) ->
	NewList = lists:filter(
			fun (Peptide) ->
				consistent_spectrum(lin_spectrum(Peptide),Spec) % remove all peptides not consistens with spectrum
			end,
			List
		),
	 % divide in two lists, one satisfying Cyclospectrum(Peptide) == Spec and one with the peptides that is
	 % consistent with spectrum Spec.
	{Acc2,NewList2} = lists:partition(
			fun(Peptide) ->
				case cyclo_spectrum(Peptide) of
					Spec ->
						true; 
					_ ->
						false
				end
			end,
			NewList 
		),
	NewAcc = Acc2 ++ Acc,
	case NewList2 of
		[] -> NewAcc;
		_ -> cyclopeptide_seq(Spec,expand_peptide_list(NewList2),NewAcc)
	end.


expand_peptide_list(PList) ->
	expand_peptide_list(PList,[]).

expand_peptide_list([H|Tail],Acc) ->
	expand_peptide_list(Tail,Acc ++ expand_peptide(H));

expand_peptide_list([],Acc) ->
	Acc.

expand_peptide(PString) ->
	[lists:append(PString,[$G]),
	lists:append(PString,[$A]),
	lists:append(PString,[$S]),
	lists:append(PString,[$P]),
	lists:append(PString,[$V]),
	lists:append(PString,[$T]),
	lists:append(PString,[$C]),
	lists:append(PString,[$I]),
	lists:append(PString,[$L]),
	lists:append(PString,[$N]),
	lists:append(PString,[$D]),
	lists:append(PString,[$K]),
	lists:append(PString,[$Q]),
	lists:append(PString,[$E]),
	lists:append(PString,[$M]),
	lists:append(PString,[$H]),
	lists:append(PString,[$F]),
	lists:append(PString,[$R]),
	lists:append(PString,[$Y]),
	lists:append(PString,[$W])].


%
% Is the (linear peptide theoretical) spectrum cosisten with (cyclical peptide experimental) 
% spectrum SpecB
%
consistent_spectrum([M|SpecATail],SpecB) ->
	OccA = find_occurences(M,SpecATail) + 1, % count this occurence
	OccB = find_occurences(M,SpecB),
	case OccA of
		S when S =< OccB ->
			consistent_spectrum(SpecATail,SpecB);
		_ ->
			false
	end;

consistent_spectrum([],_) ->
	true.




% Find number of occurences of element El in List
find_occurences(El,List) ->
	find_occurences(El,List,0).

find_occurences(El,[El|Tail],Count) ->
	find_occurences(El,Tail,Count+1);

find_occurences(El,[_H|Tail],Count) ->
	find_occurences(El,Tail,Count); 

find_occurences(_El,[],Count) ->
	Count.


%
% Test data
%

leqn() -> [0, 113, 114, 128, 129, 227, 242, 242, 257, 355, 356, 370, 371, 484].


tyrocidine_b1() ->
 [0, 97,99,113,114,128,128,147,147,163,186,227,241,242,244,260,261,262,283,291,333,340,357,388,389,390,390,405,430,430,447,485,487,503,504,518,543,544,552,575,577,584,631,632,650,651,671,672,690,691,738,745,747,770,778,779,804,818,819,835,837,875,892,892,917,932,932,933,934,965,982,989,1031,1039,1060,1061,1062,1078,1080,1081,1095,1136,1159,1175,1175,1194,1194,1208,1209,1223,1225,1322].





