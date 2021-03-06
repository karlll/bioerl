%
% 2014-01-16 karlll <karl@ninjacontrol.com>
%

%
% Pattern matching
% 

-module(pattern).
-compile([export_all]).


%% -------------------------------------------------------------------------- %%
%% Trie structure
%% -------------------------------------------------------------------------- %%

trie_new() ->
	[].

% Build a trie from a list of string Strs
% trie_strs(Strs) -> Trie
trie_strs(Strs) when is_list(Strs) ->
	lists:foldl(fun(S,T)->

			{_,T2} = trie_str(T,S),
			T2
		end,
		trie_new(),
		Strs).


trie_str(T,Str) when is_list(Str) ->
	lists:foldl(fun(C,Acc) ->

			{Next,T2} = Acc,
			trie_add(T2,C,Next)
		end,
		{1,T},Str++[$$]).

trie_add(T,C) when is_integer(C) ->
	trie_add(T,C,1).
trie_add(T,C,Current) ->
	CurrentEdges = proplists:lookup_all(Current,T),
	EdgeC = lists:keyfind(C,3,CurrentEdges),
	case EdgeC of 
		false ->
			NewNode = trie_num_nodes(T)+1,
			{NewNode,[{Current,NewNode,C}|T]};
		{Current,NodeB,C} ->
			{NodeB,T}
	end.

% Find an outgoing edge from node Pos
% named C. {Pos,NodeB,C} return if existing, false otherwise.
% If C is not found but there exists a stop edge, 'leaf' will be returned.
trie_at(T,Pos,C) ->
	Edges = proplists:lookup_all(Pos,T),
	case lists:keyfind(C,3,Edges) of
		false ->
			case lists:keyfind($$,3,Edges) of
				false ->
					false;
				L ->
					leaf
			end;
		E ->
			E

	end.

trie_num_nodes(T) -> trie_num_edges(T)+1.
trie_num_edges(T) -> length(T).

%% -------------------------------------------------------------------------- %%
%% Trie match
%% -------------------------------------------------------------------------- %%

trie_matching(String,Patterns) ->
	T = trie_strs(Patterns),
	trie_matching(T,String,1,[]).

trie_matching(_,[],_,Matches) ->
	Matches;
trie_matching(T,SubStr,Pos,Matches) ->
	NewMatches = case trie_match(T,SubStr) of
					match ->
						[Pos|Matches];
					no_match ->
						Matches
				end,
	{_,NewStr} = lists:split(1,SubStr),
	trie_matching(T,NewStr,Pos+1,NewMatches).

trie_match(T,String) ->
	trie_match(T,String,1).

trie_match(T,[C|CT],Pos) ->
	case trie_at(T,Pos,C) of
		leaf ->
			match;
		{Pos,Next,C} ->
			trie_match(T,CT,Next);
		false ->
			no_match
	end;

trie_match(T,[],Pos) ->
	case trie_at(T,Pos,$$) of  % at end of string and next pos is a leaf == match, otherwise no match
		{_,_,$$} ->
			match;
		false ->
			no_match
	end.

	
%% -------------------------------------------------------------------------- %%
%% Burrows-Wheeler algorithm
%% -------------------------------------------------------------------------- %%


% Apply Burrows-Wheeler transform to a string (Add '$' as string terminator)
bwt(String) ->
	bwt(String, length(String),[]).

bwt(_,0,Acc) ->
	S = lists:sort(Acc),
	lists:flatten(
		lists:map(fun(E) -> lists:last(E) end, S)
		);

bwt(InStr,Count,Acc) ->
	RolStr = rol(InStr),
	bwt(RolStr,Count-1,[RolStr|Acc]).

% Rotate string one step to the left
rol(String) ->
	{Prefix,Last} = lists:split(length(String)-1,String),
	Last ++ Prefix.

% Invert a BWT. Result is the first item in list, disregard terminator
inv_bwt(String) ->
	C = lists:map(fun(S) -> [S] end, String),
	inv_bwt(C,lists:sort(C),length(String)-1).

inv_bwt(_,Col,0) ->
	Col;
inv_bwt(Orig,Col,Count) ->
	NewCol = bwt_add(Orig,Col),
	inv_bwt(Orig,lists:sort(NewCol),Count-1).


bwt_add(A,B) ->
	bwt_add(A,B,[]).
bwt_add([],[],Acc) ->
	lists:reverse(Acc);
bwt_add([HA|TA],[HB|TB], Acc) ->
	bwt_add(TA,TB,[HA++HB|Acc]).

%% Matching using Burrows-Wheeler

test_bwt_match() ->
	 S = "TTCGCTTCACGTTCATCATTCGCATGTCTAGAGGATTGATGAGGATATCCACAGCGTCTACCTCAACGGATCTAGATAGGTAATGGCCGTCTCTTCTAGCATCCGCCTGTTGTTGAACCCGCATGGAGCACACGAAAACATAAAGTCCGTTTTAAGCTCGATGGCCCATCCCTGCCAAATGACGCAAACCGTATAACTAGACTAATGGTGTAGTATCGGGCGTTTCAAAAAGTACTCTGAGAGAAAAAATGAGTTTGCTATGCACTCCTAAGTATCAATCAGTTCTGGTGAATCATACCGTCGACCAGTAAGTCATGTTAGTAGTTAGTTCCAAAAAA$TACTAGGGTGCCTTATATGGTGGCTCAGTCGCTTGCCTTAGAGCGCCCCAATGACCATCTGCTTCGCTCGTAGTCTCGCCTCCATAGTGTAGGCCCGGTCGCGAGCTGCAGTTCCAGATGGTTACTATGGTGGCCTAGACCTACCGCCGATATATGTCCGATGCGATAAGTTAAGCATGCGATGGACGCGGGTGCAAAGCCTAGAGCGAGCCTTGAACTGGACTACTTTCCACAAGCAGTCCGCACAAGTTTACACGTCGTGGGCCTTATCAATTTGTTTTTAAGGAGGAGATTATGTGACAACGGACTGTATAGCATGGTGCTAACTCCGTTGCCTGGGCCCCCACCCTCAACCAGTACTAGGGCGGGGAGGAATGTCATAGTTCGGTTGCGTGTGATGTATTCAGGTCACCATGCACATCTGAAGTATACCGACAGAATAGCCCACGGGAAGAAGGACCCCGACTAGGGGAAGCTTGTCCCGACCAATCACTCTAGAGAGGGGAGATGTCTTCGCAACCCATACACATCAAAGTTATGATCTAGGCCCAATGTTCAAAAATAATACAGGACGTGCTCTATGTTCGGCCTAGAAGAATAATGTCATCGACAGTGGGTGCTCTGGCAGTAACTAGCCTGTGTGGTTCCATTAATAGAACATC",
     P = ["CGACCGGATG","ATACATAAAA","CATTGCTTTC","TCGGGACCGA","AGCCCCTAAA","GCCAGTAAGA","CACTTGAGGC","ATTACAACGT","ACTCAAGTTG","TTATTTCACT","TTATTCGGTC","TAATGCCCCA","TTTGATAGTA","AATTATACGT","TACCCTCTAC","GTGTGCACTT","AGGGTACGAC","CTGGACTGGA","ACAGCCTATC","AGGCCGGCGC","TTTTCCCATA","ATCTCTATCT","TCGCAATCTT","GCCAGAAGAC","CCAAACGTCA","CCTCGATTCT","AAACCAGGGC","ATAATTTCTA","TTATATTCAA","GGTCAGACAT","GGAGACTTCG","CAATCTATAT","CTTCCACGAC","GGTGTAGCTC","CCTAATCTCG","CGGCGGTTGG","TGATAACGGA","CAGAGATTAA","AAATTTGTTC","TTCTCAACAC","ACCTCAAAGG","ACACGAGCAT","ACCACCTCGA","GGCTCGATGA","TCAAAATCTA","AAGTAGGGTA","AATGTTCATT","CTTCTAAAAG","CATGCGTAAG","TCACGACTAT","CTGACAACTC","GATTAAGAGG","GTGGGCGGGA","TAACTGTTCG","GCTTAGGGTG","ATATCCTCTG","TCAAAGGACA","ATGCTAACTG","TGTTAGGTCG","CTATGATAGT","CTAGGCGCAG","GGACTGGAGT","AGTAAGATAT","AACTACGAGA","TACCGCCCCT","AAAGCTATAA","CTAAAGAAAG","TGTCGTACTG","GGCTTTGTGA","TTGGAGCATC","AAAACTCTAA","AATATATAGG","AGTTCCACTG","TAGAGTCCAG","GACACTATTC","GGGAGTTGCC","TCCAAAGAGA","TTTCGAGGTT","AATACTGAAG","ACTGTAGCAT","GGACTCACTT","AGTGAGCACT","GTACCAACAT","GAATGTAGCA","ACACATGCGT","ATTCTCTTCC","AAGACACTAT","GGGCGTTGTA","TCTGGTATTA","CTAAGCTAAT","CACAGAATGC","TGGCTTGGCC","GACTAATTCA","GGGCGCTATA","GTCCCTCTAA","GGGATCGATA","TACCAAAAAC","TGTGGGAATG","TACTGAAGAC","GATCCCTGAC","TAAATGTAAT","CCCTGCTATA","TCATAGTTAT","TCCCTGCTTA","TGCGCTGCCG","AGCGATCACC","TTGAGGATGG","CCTAAATGCA","TTGGACTATC","CCGAGTGAAC","CGCGTTGCGA","CCTTCCATGA","TCAAAAGCAG","TCTATCAGGC","GGGAATAGGA","TACGCAAGCA","TGGGCGTATC","CCAAACGATA","CGCCGGCCGA","AGATCATGCT","TACCCCGTTT","ACCGCACTAC","CGGTTCAAAA","CCTCGATACC","GGCATACGAG","ATAGCCGGTG","TCCGATGCGG","TCGGGTCATG","ACAGTCGAAG","AGCCGGTGTT","AATCAAAACC","TTGTCTTTTT","GTGAAGAAGA","AACGAGGATA","AGCCTCCTTT","CCTTCTTAGC","ACTTCCGCTA","ATGAGGACGG","TGCTAAAGAA","AAACCTCCAC","ACAGCATAAA","CTTCATACAT","AGGATCGGAG","AAAAGCCCTG","ACATGGGCTA","ATTCCCCCAG","ACCTAAGATT","CAATGCTCAT","TATTCGAGTT","AGGCTCCTCG","AATACACCCT","CTGCTATACC","CCTTCTAAAA","ATTCTGGTAT","ACATTTGTTA","CATTGTGATC","AACCACTTTT","TCGATTCGGT","CTCCGGTTCA","GAAAGTAGCC","GGCGGTTGGT","CGGTGTAGGA","CTACGAGCTA","AGTCGTCCGG","TACAGAGTAT","CTTTCGGTGA","TCACCACCTC","GCGCATCAAT","GATGCTGGAC","TCATAATCAT","AGCGCCGACT","TGTACATGGG","TCAAGTTGGA","TGTGAGGAGA","GCGCCCAGGG","GGAGCTTAAG","GGGTCCAAAC","CTGGGGCGCT","GAGGAGGTGT","CTGACAATGG","GCCCAGGAAT","AGATTCAGAT","CACAAATAAA","AATCTCAGTA","AAGTGTAGTT","CCGTGAAATT","TGGGGCGCTA","ACCAATTGGC","GTTATTTCAC","GATAACCCAA","TGATCCCTGA","TTTTACCTCG","CTTTTGAGAT","GTGAAAGGAA","GTGAAAATGT","GAATGTCTCG","AATAATGTCC","ACTGTTCGCG","GTGTACGTAG","CCCGGTCGGG","TGGGGAGCGC","TAAACCTGTA","TATCGTATAA","ATAGTTGAAA","AATGCTTATC","AATGCTTATA","GGAGTCTCTT","AGCAATGACA","AAGAAGGGTG","AACTCGGGTC","AATTTCTCGC","GGGGCATCAT","TTAAAACTTC","TGGAGTCTCT","GAGTTTTCGC","CCCTGCTTAA","TTCGCGTGCT","TCCGATGGTC","AAAAGGAGCA","CTGTGGATAC","TAGTATGTCT","AGCAACAAAT","TCCCTAGCCC","TTAATAGGGA","GGATCGGGTA","GGCTCCACTG","AGTGATGTTA","CACAGGATGC","AACAATCCGG","ACTCGTCTGT","CAATCACGAC","CGGCACCAAG","CGCGGTTTTC","CCAGCTCACA","TCATTGCAAG","CCATTATGAC","ATTTTTCGAG","ACAATACTGA","AGGAGGTGTA","ACCTCGTTCA","TTTGTTCAGT","TCCGCCTCAC","GAGACCAGAC","TCCAGTGTAG","CAAAAGAGTG","AGGTAAAAGG","TGGTTCACCC","TGAAGAAGAC","TTAGCAATGC","CCATAGGTGC","AGTGTCAGTG","ATTGATTTGG","GCTAACTGTT","GTTCGGCATG","CACCACGCTC","TGTGCTTTTC","CACCCATGGG","ATGGTTTAAC","GCAGCGTCAC","ACAGGGACAC","ATCGGCCCTA","TGTATAGCAA","CATGCCCTGT","AGCAACATTG","GGTATATGCC","TATGCCTGCC","GTAGTTAGTG","TCCCCCTATA","TGTGATTTTA","GCTCACTGCC","TGGCCACGTC","GTAGCAACAA","CGCAAACGAG","TGACACTGTA","CATAAAGAAT","CCGGGCCACT","CGTATTGAGA","AAGTCCCCAA","CTCACAGGAT","CTCTGGGCTG","GGAGACTACG","ACTGGTATAT","GAGTTCAGGC","TTCTCAAGGT","GAGATTAAGA","TGCCCGGTAA","AGCTTCATGG","TCGGTGGGCT","GGTATCGGTA","AGTAGCCCCT","ATCTGTAGAC","GATTCTTCAT","CAAATACACC","GAGTAACAAT","CAAGATCGGA","ACTTGCTAAA","CCTGCTTAAT","CCAGCCCAAA","AGTAAATCGA","CGTTTGGGCG","GAATGGCGTA","CAAAAGATAC","ATGCCCGCGT","TTAGTGACGG","TCCATATCTA","GGTATTAATC","CACGAGTGAA","GGCCTCCTAC","CTCCAGTTCC","TGGACCAAAG","TCTACAAATA","TAGGTGAGTT","GCTGACTCGC","GACTAGTTCA","GAGGAAGCGC","ACATCTACAA","GTTCCTACAG","AGGCCCTCGC","GCCGGTATGC","AAAGAGAAGT","TATCAGGCCC","CAATACTGAA","AAGACTTGCT","AATCTGGGGC","ATAGGGAATA","GGCGGTGGGC","GTTTCCTTTC","CCCACTTACG","TATCGTACAA","TAGCAAAGAG","CGGTTACCCC","ACTTATCCAA","GCGGTCCCTC","CCGGTGTTAT","GGTCCTTGTC","GCCGACTGGT","GAGGCCTATA","AAGAGGAGGT","CCCAGGGGAT","GCCCCTAAAT","TCCCTGACAA","CTCGGTGAAA","CCCCTTCATA","TCTAGGCTCC","TGGGAATGTA","AATGTATCTT","GGGGCTCACT","CACTACGGGT","GCCATCAAGC","GCATAAAGAA","CTGCTTAATA","AAAGGGGGTA","AAATACACCC","AAATACACCA","CCCTTCCGGG","CGCAACGCTT","AACGTGGTAG","TGAGGACGGT","ATAAGTAGGG","TCTCAACTTA","GTATTTGTGC","GGTGTTATTT","ATGACAGACT","TTTGGCTGCA","CTATGGCACA","ACGACATCAT","ACTAAATTTG","GTGCGCGTAA","AAAGCGCCCG","GTGACCAATG","AAAGACATAC","GTATGGACAA","GTGCTTTTCC","ATGAGACTCT","CTACCCGCTT","TATCATACTG","GAAATACACC","ATGCGACTGC","GGGCGCCACG","CACCAACGGT","GATACATGGA","TATATAAGCG","ACGCAGACTC","GGTGGGTATT","GACGCCGGCC","TTTTCCGTAG","GCTAACCTTG","ATGTCCCTGC","ATAGGACCTT","CCAGTACCAT","CAATCTGCAT","CCCGGTTACC","CAACAGATCC","GGATAACGTG","TGCCCATATG","TGAGCTGCGA","TGGCGTAAGT","GGAACGTGGT","AGATCCTTCA","AGCTTAAGAG","ACCTTCTAAA","AAGCTCGAAT","GCACTTTAGT","AGATCTGGTA","AATAAGTAGG","CATATCCGAA","CATTGAATAA","GAAGAAGACT","TAAGCGATCA","TCAACGCCTC","ACAAATGCGG","GTGGACCCCT","GACTCAAGTT","CTCAACACCC","CACACATGCT","TCGTACAATA","TCAATACCTG","GGGTCATACT","CAGTATAAGT","CGATAACCCA","CACCTCAAAG","TGCTCCGGTT","CAACGTTGGG","TAGGATCCGA","CTTACGTGCA","TTATTTGCAG","GTGGTACTGG","CTGTACTCGG","CGCCGCGGCC","TTCATGTTTC","TCCCGCTACG","TTAAAGTCCC","ACAGGACGAG","TCGCGTGCTC","AACGTCAGCC","AAACAGCGAA","TTCCAACTAG","TTCGTCGCAA","GTTTATGATT","TTGCGGTTGT","TCACCTCAAA","CACTGGATGG","TCATAACTTT","AAGATCTTTA","GTATTTGATG","AATAGGGAAT","AGAAAGTAGC","TTGCGACCGT","GCCGCGCTAG","TGTGGTTTGT","TTTATATATG","AGCTCGAGCG","ATCCGTCATT","CGAAGTTTCC","TTCTGGACCG","CGTCGTCCTG","AGGTAAAAAC","TTATAGACTA","ACACTTTAGG","GTGGTAGTCA","TCGAGACTAA","GATAAAGAAG","TCTCAACACC","CTACCAGTAA","CGGAAACATT","TTGTGTAGGG","TTCCGTAGCC","ACGCGGGAAC","TGGGTATTTG","TCGCTGCGCT","CTCGGTGATG","ACACTATTCT","ACCGTTGCGA","CAGCCGAAAT","ATTAATGCAT","ATCAGCGTTG","TCAATGCTCA","CGTTATTCAC","AGTGAACACA","GCACTACAAA","GTGCTCGCCG","TGCTAGCTCT","CTTAGTAACG","TAGTCTTTCG","ATCCGACCGG","AATCATCGAT","AAGATCGCGC","AAGATGACAA","TATGAATTTA","CAGAACATAT","AACATCATAA","CGGGTTGTCC","CTTCCGGGAC","AGGTGTAGCT","ATGACACGTA","CCGTCATTTG","TGGCATGTAT","AAATCTATAA","ATCTTTTAGA","GATTGCCGGT","GCAAACCGAG","GGTGCTGGTC","CTAAATTTGT","CGGCCGGGCA","CTCGCAATCT","AGACAGTAGG","CAGGCCCGCT","AGCATGGCGA","GTGCTACTTA","GTGCTCCGGT","ACGTGGTAGT","GCGAATGGCG","TCATACATAG","GACGATCCCT","CTAGAAACTC","TTGTTAACTA","AAGTAGCCCC","CGCCAGCACG","CGTACAATAC","CGACTTGTAG","ATACTGAAGA","AGTTTTTGTA","ATGCGAAAGC","CCTGCAGGAC","TCCAAAGGAC","AGGGAATAGC","AGGGAATAGG","TATCTTTCTC","CACGTATTGG","TGCTTAATAG","ACTAACATGT","GGGTACTGCC","CCGGAGTCAC","TAACTAAATT","GAGTCTAGGC","CGCTACCGAG","AGATGACTCA","GCCCGGTTAC","TGCTCATAAA","TCCTTAATGA","CCCGGTAGGT","TGCACGTGCC","CCAAAGAGAA","GTTGGATATC","ATGGAATCTT","CTAAAAGATG","GGACGCCTAG","GGAATATGAG","ACCTCAGATC","ATGCTCCGAC","ACCTCAATCA","CGTAGCCGCA","CAGATAAGGT","CATTTACACC","GTCACGATCC","TTACCTATTG","TCGGCCGAAT","AGAAGATCGC","TAGAATATTC","TAAGATATCA","ACTACGAGAG","TGATCGTTCT","GATTGACGGA","AGGAACTTAC","GGAATGAGAG","CCGCACTAGT","TATGTTACGA","TTAACTAAAT","AATATCGTTA","CCGGTAAGCC","GGACACTACG","CGCAGCAGGC","CCTACAGGGA","GAGTCTCTTA","AACGCCTAAA","TCAGTAAATA","CTACGAAGGC","GTAGCTCAGA","CCCTGACCGA","TTTTGGAGCA","CGCCGATATT","GGAGCATTTC","TTATCCATAA","CACGGACCTC","CTTGCTAAAG","CTAGGCCATG","CACCCTGCTA","TTTTAAGGGA","GTCGAAGGCG","AATGTTCTCA","TCATGCCTCG","GGGGCGCTAT","GACAATAATG","CCTAAATCCA","AAGAATTTCT","GTCTCTTAAA","TGTGATCCGA","AGACGCGGGC","TGGTTGGGAT","TATTCCACTA","TTCGTAATTA","ACTCGCCCGA","GCGCAGCAGG","CAGACGCAGA","TTGCGAATGG","CAGGATGCTG","CGCCCGGTTA","TGCCCACCCA","GTATAGTTTC","ATGATGAGGT","CCATCTCCTG","GCTTCTACCA","CTCTTCCCGA","AAGGCGATGG","TCTCGGACCG","TGCTTGTTGT","GCCGATAACC","AGTCAATCAC","ATTTGTTCAG","CCACAATTTG","GATCATGCAA","GGCTAGAGGG","TATACCACAG","TATGCCCGCG","CAAAGGACAA","GCTCATCACA","AATAGGACCT","AGAGCGCATA","ACTACCGCAC","GGGAGATAGC","GCCCGCGTTG","AATTCTAGGT","CCAAAGCCGT","TCTATTCCCT","GACACGTAGT","ATGTAGCAAC","GACTTATCCC","TAACCAAAGC","TTAACATCGG","ATACCTTGTG","ACCCTGCTAT","AAACTCCGGC","TACAGGGATA","TTACGTTATG","GGATTTCCGT","ATTGGCAGTT","CGAATTCTTG","ATTAATCTGG","TCAGCCTCAG","TGTGAAGACG","GTGTTACAGG","ACCAGCTCAC","CTTAGTCACT","TATCGGCGGT","ACAACCCCAT","GCTACGGCAT","CGTCCGCGGT","GTGTAGCTCA","TACGTGGCGA","AGGTTGGGAA","AGACACTATT","TTGGGGCTCA","AATGTGATCC","ATGTGTCTAT","AATCTCTTAA","CTATAAAGGA","GTAAGCCTTG","ATACAATTCA","TTCTAAAAGA","GCGGGGCCTG","CCTCGGTGAA","AGTGTTAAAG","TCAACACCCT","TTCCTACAGG","TGGAGCATCC","CGGGCCCACA","CGGCGCCTGG","CGATCAGTAG","CCAGATACCG","TTCCATTTGT","ATGCCTCGTA","TGTAGCTCAG","TAATTGACAG","CACTGGTGCC","CAATAATGTC","TGGACCCCGA","TTTTCCCCTA","ACGGACTACG","CAAGTTGGAT","TCCACATAAG","TTGCCGCATT","ATCAGGCCCG","GGAATAGGAC","TTTCCGTAGC","ATGGAAATTA","TCGAAGGCGG","TCGCAGTTAC","ACCGAGGATA","GAATGCTCCG","TATGTATAGG","ATCAAGAACG","ATCTCAGTAA","GCTAGCGATT","TGTCCCTGCT","CATGCACCAG","TCGCTTCTCA","TCAGAACCAG","GTGATCCGAC","AAGTACATTT","TGTGTACTTG","AAACCTACTT","GTACTGGGGC","CCTTCCGGGA","TCCCGGTTGG","CCCCAATGGT","CGGTGAAAAT","TGTGTCTATC","CAGCATAAAG","TGCGCACTCC","AATCATGTAG","TCCCTCAGCG","GAGATTTAGG","TAAAGAGCGG","ACTGAAGACA","TGAAATATTA","GGCAGTTGCT","GCTGGGCTAG","GCTCGTCAAA","AGGGTCTCCC","CCTCGGCGCC","CTGGTACGTA","GCTCAAAATC","GAAGCGCCGA","GGTCATAAGA","TGCATGGTAG","ATAATACCAA","ATCGATAACA","CCCAATCTAT","CTATAATCAT","TAAGCGAGTG","TCATTGATTT","TACCTTCTAA","TTATTCTGGT","GTAGGGTACG","AACACATAGA","CCGGGACGCC","ACCAATGGAC","AAAACGGTTA","ATCTATAATC","CGCCATCAAC","TGACTGGTAT","GGGCGTATCG","ATCTTCACGA","CGCCAAATCT","GCACATCTAC","ATAATCATCG","GTTATTCTGG","ACGGTGCTCG","TCTGGGGAGC","ACTTACAACG","AACCCTCCCT","TTTGGGCGTA","GCTTAATAGG","GAATGGCCCA","CCTACGAGCT","GATCGATATG","TATGAGGACG","TCTAAATGAC","TCATCGATTC","CCTAGGGACT","GGTCTGGCCG","TAAGTAGGGT","CTGTAATGAT","TACAAATACA","TCCGTAGCCG","CTTCTGGGCT","GCACCAGCTC","CTATCGCTAG","GATTCTTAGA","GAAGGCGGTG","AGTCCTATCT","ACCAGCAACC","CGTAAGTAAA","CCATAACACG","CCTCCATATA","CTCTTAAAGT","CTACAGGCTT","GAAGATCGCG","TGTACAATCA","GCGATCACCT","CGAGGCGCGC","CCGTCTGGAT","CGTCGTGTTG","TCCCAGTTTA","TTAGAGTCGG","ACGGATCCCC","TATTTGATGC","CGTCCCAACA","CGGTGTTATT","CCGCATTAGA","GCGAGGTTGA","CCGATAACCC","AGGAGCCCGT","TAATGTGATC","CGGCCGATAA","CGCAGGCAAA","ACAGTGCTCC","GACAAACCTG","GACCAGGCAG","ATCATTGATT","CAGGGGTGGG","TCTCGTGGCG","ATTTCACTAG","ACCTCGAGGA","GCATTTTCCG","GGCGGTAGCG","TCGGTGAAAA","GTCGGACGTG","ACTCGCATTT","CCGCGTTGCG","CTACGACGTC","AATTTATGTG","TATATAAACG","ACTATTCTCA","GTGTGGTAGA","GGCCTATGAG","CCCTCTAATG","CTTCTAATGT","CTGTACCGGG","CCCATTTTGG","TATTGCCCAC","GTTACCCCGT","GACACATGTT","TCAGGAGATC","ACTTCTAGGA","CACGTCCTTA","CTGAATACCA","CGTAACATTT","CCTCTAATGT","GTGGGTATTT","GGAGAAGATT","TCAGAAGATC","GCGCCGACTG","TTTGTTAACT","CCTTAATAAA","TCTCTAGCTA","ACCACAATTC","CACGTCGGTT","TGCTATACCA","ACAGAAGAAA","GCCTCTTAAG","CAACTCGGGT","TCTATATACG","CTCTTACTAA","GCTTATAGCG","GAGTTGTAAG","ATTTTGGTCC","CAATGCCAAC","GGGCATATTG","GCGCTGTTAA","CTTCGTTGAC","GTTTCGTCTT","GGTTACCCCG","ATCGGATCGT","GAGGTGTAGC","AGTTAGACTC","AGGCAGAGAA","ACTACTAACA","AGCTGGCCTG","TGCTCGCCGC","GAAAATTGAC","TACAACGCGA","GGACCTTTAT","AAGAAGACTT","CCGGTTACCC","GCTCCGGTTC","CCTGGGCGTG","CTGCAATTGG","GGCCGATAAC","TCCAAACCAC","TTAACAGCTA","GAAAGTTGCC","AGACTTGCTA","CACAGCATAA","CCTTCCCGCG","TCAAACTTAA","GCGGTTGGTT","CACGCAGGTA","TTAGTATGGC","TATTATACAC","GGAGGATCCA","AAACCACTTT","GATCTGATAA","GGAGCATCCA","TACGAAGGAC","CGCTATAGCC","AAGTTCCTAC","TGCCGTGAAA","TTGCTTGTTG","GGGTCATGCC","ATCAGGGGCA","AGATAAGGAT","TGTTCTCAAT","CCGCTCTGGG","TCGGAAATTA","CGCCGACATA","GCTGCACAGT","GGTCTGGTCC","GCCTCGGTGA","ATATGTGGGA","CTCCAGTGGC","ACCCATATAT","CCATGGTCCA","TGTTCCTTGT","CAAGAAGGAG","GGGGACCAGA","CCCCTGAACT","GTCGGCGACA","TTGATCCCTG","TTGACAGGTC","ACTAAACTTC","AGGAACGTGG","ACATCCGATT","TTCGAATCAG","GTACGACAGG","AGTTATTCTG","AGATCCCACC","ATAGCCGTCC","AACAAATGCG","CAGGCTACCC","TGCGTTAATT","ATACGGGCCT","TGCTAGACTG","TCTTCGACGA","TTCTAGATCG","ACTTCTCTCT","TTCTCAATGT","GTCCCGGGGC","TCTTAAAGTC","CCGTACGATT","ATGTCGACAA","CGAGGGTGCT","ATGGATATTT","GATCCAATAG","AACTCGCCCG","TCGGCGAAGA","GTGGAATAGG","TATACTGGCG","CCGCGGTCCC","CGTTGCGAAT","CTATGGTTAA","CCCCGGCGCA","GGTTGTAGTA","TAGCTGGACA","GTCTTGCCAA","AATGGCGTAA","GATATGTGGG","TGGCTGCACA","TGGAGCCTTA","GCAATGGATC","GACATTCCGG","ACCCGCACAG","ACGAGCTACG","GGACGCCGGC","GAACAATCTC","TATGTGGGAA","GGGAACTCTA","TAGCAACAAA","ATCACGTCAA","ACTGTTTGTT","AATGTCCCTG","GGGATACCTT","ACTGCCGCGC","AGACGCAGCT","GATCGCTTTT","TGCGATGGTT","GATCTAAACC","CTTAGTACAG","TTAAGCGGGG","GGGTATCGGT","CTAATGTGAT","CTCGAGGAAG","TTATCACTAG","ATACCTCCCT","ACTGCGTTAC","CTTAGACAAC","CAAATGCGGT","ACCAGAGATT","GTATGAGGAC","GGACATCACG","CAGTAAATAA","CGTCAGGCTC","CACACCTTTT","CATCATAATG","CCCTAAATGC","TGACCCCGAC","GGGATCCCGG","GCGTCGTTTC","CGAATTGATA","CGGTGCTTGA","ACACCCTGCT","ACTGATTCTG","CTCAAAGGAC","AAATTACTAC","AGAACAAAGA","ATTGTCTATA","ACCGGATGTG","CAGGAACGTG","CGAGGAATCC","CTATCAGGCC","AATAGCTTCC","TTCTGATAGT","TACGCAGCAA","ACACTCATTG","CGATATGTGG","TAGTCAATCA","TAGTACACAA","TCTCCGATGG","CTCGGATCTG","AGCTTCGAGA","TGCTTCCGTA","CAGCAGGCAT","TGGCCCGTTC","GAAAGCCAAC","TTCACTAGAG","ATGCTGGACT","CTGAAACCTA","ATCAAGGGAA","AAAGCCAGGG","TAAATGTTCT","GGGACGCCGG","TTGCATCTCC","CGCGCAGCAG","TAAAAGATGA","TGATGTGGAC","TTTGGGGGCT","TTATGTTGTA","CCTCACACCC","GATCTTGGAT","GGAACCGCGC","ACAAATACAC","AGGGCAGAGG","CATGCTTGCG","GACATGTCTC","TGCGTACCGG","CTTCTACTTC","GTGTGTTAGG","CTTAATAGGG","GGAAGCAACG","AACAATCTCA","GTGACAAAGA","CCAAACCGCA","GGCCCGCCCG","CGCGTCGATC","CGGGACGCCG","TTGGGATCGA","CATTTTCCGT","ATGTGAAGAA","GTAATATAAT","GCGACCGTTG","GCCGTCCGCG","TGAATACCAC","TTGGGGCCGT","TCCCCAATCT","ACATCAGTCA","CGGTGCTCGC","AGTAGGGTAC","ATGTCCCGGG","AGACCGCCAC","TGGTATATGC","TCCATGACAT","GGAAAAGATC","AACATGCATC","TGTACCAGGC","AAGTGCACAG","CAGTAAGATA","AGTTCATCCC","GTAGAGCCAG","CCTCGTCGCA","TATAGCCGAT","ACACGGTGAT","CACCCTTCCG","GTCTTGCATC","ACTCGAGTCG","TTTATCATTG","GGTCCAAACC","AATGTCCCGG","CCTGTCAGTA","TCACGCCAAG","TGATTGGCGT","GACTACAGCC","CGGGTATCGG","AAAGTAGCCC","GACCGCTCGT","CCCCACGCAA","GGGACACTAC","AGCAGGCATG","TAGCCCTAGC","AGATTAAGAG","GGATAAAATA","GTAGTTTTTC","CCATGACGCA","ACTTTCTAAG","CATCTACAAA","GTGGTAGGGA","AGGGGTCCGC","GTCCAAACCA","AATCTTGCTT","CTTTTCCGAT","CAAACGTATT","CCGGTTCAAA","ACGTCGGATA","GGACCAGAGA","TAGCCCCATA","TGGATGTGGA","ATCGTACAAT","GCTAGGTTGC","AAGGCGGTGC","ATTGGCTCCA","ACGGGTCCAA","TTCTCGCAAT","AGGCAGCCTT","CCTTGTTCAG","CGTGGTAGTC","TTGGTTGGGA","CAAAGCCCAC","CGAAACCACA","TGCGTCGCGA","GCTACGGGTC","AGGACCTTTA","AGGACCTTTC","GCAGAACAAT","GTGGGAATGT","AACGGTAAGT","AAGAGATTGT","ATTTGATGCT","GGGGCCTCCT","AATAACTAGT","GGGAGGACAC","TCTAAAAGAT","TTTCCATCCT","ACACTACGGG","TCAGCCTAGT","TTTCGGGAAC","GAGTCCAGTA","TCGCCTAGAG","CTACAGGGAT","CTCTCCTGAT","ACTTACAGCA","TGAGATCACC","TATATGATCG","GGGAAATTTC","TTGAAATGAA","ATGGGGACTT","TGTTAGCATC","TGACGTCAAT","GCACAAATTG","CCAAAAAGCA","ACTCCGACCC","TCTGGCGCTC","CCAGCGGCGC","GTGTATAAGC","AGATCCAGTT","TCCTGCGAAC","CGGTTGGTTG","CCTCAAAGGA","GTCCCTCTTT","ATGCCTGAAT","TCGTCAAACG","AATATAATCA","GTTCAGGTTC","ATTATGGTGG","ATACACCCTT","AAGATGCCCG","ATGACTCAAG","ACCGGACGAT","CCGCAGAACA","TTGATTGTTT","GACAGGGACA","GAGGACGGTG","ACCCGGCTTC","GCTATCGCGT","ATATTAATCT","GTATCGTACA","CATATAAATA","TCTCTTAAAG","CTTGATCCCT","ATTCTTCATA","ATTTGGCCTA","CAATCTTGCT","GGCATGCACC","GGAAGGGTGT","GTCGTGACGT","CCGGATAACG","CTTGCGGAAA","ATTAAGAGGA","CTACAAATAC","ACAATTTTAG","GATGCGGTAT","ATTCGTGAAT","TCTCAATGTG","CCAAGTCTTA","ACTGCACATG","GAATTCTCCA","GTCATGCCTC","CTCTATGGGC","AGTCCATTTA","GTCTAAGCGT","AAAGGACAAT","GTGTCAGGAC","CTCGCCACCT","TTAAAAGACT","CTACTGTCCA","CTGGTATTAA","GTGACTTTAA","CGTGGTAGGT","TTGAGTGTGG","TAAATTTGTT","TGACTCAAGT","GTTCTCATCT","TATGCGTCGT","GTCAGATACC","ACATCCCATT","CACTGCCAGT","GGGACCAGAG","TCGCCGCAGA","GTTCCTAGCT","TTTGCATCAC","CTGGCTGTCA","ACTAGAGTCC","ATACCACAGC","CTCGCCGCAG","TCTATACTTG","AGAAGTACAT","ATAAAGAATT","GCAACAAATG","GGGGTGTTTC","GGGAGTAATT","CTCCCTTTGT","ATAATGGACA","TAGACGCACG","ATCCCTGACA","TTAAAAGAAC","TCTTGGGCCT","ACTATCGTCC","CATGTTACCC","ATGCCTCGGT","TCGGCACTCC","CAACACCCTG","CAAAGCCCTG","GTTGCGACGC","GATAAACTCC","GTAGTCAATC","ATGCACATCT","TCTGCTTTGA","CAGGCATGCA","ATCCAAAGAG","GGGGGGTCGC","ATTCGATTGC","ATAAATGTTC","TGGTTTCTTA","TGGACTCAGT","TATTGCGACG","TGTTGGAGTC","TCGTAGGGTC","TACATTTGTT","CCGCAGGGGG","TCTACCATTT","AGTCTCTTAA","ATAAGCGATC","TGAGTTCTAT","CGCATTCACC","CGGGCTCATA","TTGTGCTTTT","CGGCATTAGG","CTTTTTGGAG","GGTAATGTGC","TTCACCTCAG","ATAGATTAGC","ATGTGATCCG","TGTTCAGACG","TCAAGGGGGG","TGCGAATGGC","GGGCTTCCGC","ACTCATGTAA","CGAGCGTCCC","GACAATGGGA","GTATATGCCT","AGCACAGTAC","CGGATATGAG","CGGGGTCCGA","ACCTGTTCAA","CATGCAGACA","GATGCCTACT","AGCATAAAGA","TATCAGAAGA","ATATTATCTT","GGTTCAAAAT","GCCGCAGAAC","CCCTGCGCAC","CCGGCGTGTC","ACTCGGGTCA","CCCGCGCTAG","ACAGGTGACG","ACAGGTGTGG","GAGAAGTACA","TCACTAGAGT","TCATCGCCTA","AGAGGACAAC","TAGGACCTTT","CCACAGGTAC","TCATAAATGT","CCCCAATCTA","AAGAAGGTCT","TAGTTTTGAC","AGCCAGCTGC","CATCTGTGAC","TCCGCGGTCC","GCATCCAAAG","CCTGACTTTG","TGATTTGGCT","ATAACGCTAC","GTCCCCAATC","CGCCCTGTAA","AGCGCAGTCC","TGCTTTTCCG","TTCCGGGACG","AAGGCGGGGT","AATACCGACT","CACATCTACA","GTATTAATCT","CTCCTAAAAC","GTCCTTTCTC","AAGGGCTTGG","GCGTGCTCCG","GGAGACAGGG","CATGACGACA","TGATCATATC","CTGTAGATTC","CATGGCTAAA","CCAAACCACT","GAGGCGACAG","GGTAGGTTAG","AGGACGGTGC","AGAATTCTTT","GCTCAGAAGA","ATCACCTCAA","AGCAACGGAC","CAATGTGAAG","TCTAGACCCT","GGTGCTTGAT","AGGCTCTGAA","CCACTCCTTG","CACTAACTTG","CAATCACCAC","GGGAAACGCC","TCGCTTTAGG","GCGCTATAGC","AACTGTTCGC","GAACAACCGC","TAGCCGGTGT","TATGTACCAT","GCCAGCACTT","CTCTCTGAAT","GCGCGCAGCA","GTAGCCGCAG","GCGGTGGGTA","GGGGTACATA","AAAGAAAGTA","CACGACCAAG","AAAAGATGAC","TAAGGTTCTG","CTTCATAATC","TCCCTCTAAT","AGGCATGCAC","AAGGGAGCGC","CAAGCGTAGA","CTATATAAGC","GTGCTGGGCG","AGACAGGATA","AGCCCTGTTC","CTTCCTTGCC","CCATGGTCGA","GGTCGCTGTG","TAGCCGTCCG","TAGCTCAGAA","TACACCCTTC","GGATGTGTCT","GGATCAGACG","GCCTATCGCT","AGTACTTCTC","GCTAAAGAAA","TTGCATTCGT","CAAAGAGAAG","GCAGCGCTTT","CCCGCAGGGA","GCTCGAAAAT","CAGGGACACT","GGTGAAAATG","ACTCTACCCG","CTGGACGTCG","GAGCTACGGG","GCGGTATGAG","ATGATCAGCG","TTTTCCGATG","GAGGGGGGTT","ACATTCGTCT","TGGGGTAATG","AAAGTCCCCA","ACGGGTATCG","TAGTCCCCGG","CTGGTCCATA","ATGGCGTAAG","AGTTAAGGTC","CCACCTCGAG","ATACCTTCTA","CGAGGAAGCG","GAGTTTACAC","ATCATCGATT","GGTTACTGAA","TCATCACGAG","TGATATCCAC","AAGAAAGTAG","AGACTCCATG","GTACAATACT","TACCCAAGGT","TGCGGATACG","CTGACCACTG","GGACGATGAT","AGTCGTTTGC","CTTTTTCAGC","GGGACCCGGG","CACAGTCGAA","GGTGCTCGCC","CGAAAAGATA","GGGTGCTCCA","ACGTGCGACA","TAATAGGGAA","GACTGGGTTG","ACGAAACCAT","GTACATTTGT","GGACGCACAC","AGCGGGAATG","ATCGATTCTT","ATTGCATGAT","AAAGAATTTC","TAGGATAGCA","ATCTCGCTCA","GATGCAGGAA","TTCCGATGCG","TGTCCGAAAT","TGAATATTCC","AGCACGTTCT","CGGAGCATGT","CTGTCGAAGC","ATACCACTCG","TGATGCTAAC","GTGCTTTGCT","TGCACAGTCG","GGGTACGACA","GCACTGGGCC","GGTATGAGGA","CCTTAGATAC","GCCGAGCTCA","CCTTTATCAT","TAAATTCGTC","CCTCGAGGAA","CCCGCGTTGC","GAAAATGTCC","GGACTAGCAC","GCTACAACTG","GTCCAGTAAG","TCTCCAAGCT","ACGACCCTTC","CCCGCCCGGT","GCACAGGTTA","AGAAGACTTG","ACAACACGTG","AGAGGAGGTG","CTATACCACA","TTCGAAGTGT","CAATGGGATC","TTTATGGGAG","TATAAATTTA","TAGTACCTAC","GAACTTCTCT","CATGGGTACA","GCAGTACACG","GTAGGGAGAT","TACAATACTG","CCTCTCTCGC","ACATAGCCGG","ACAATCTCAG","TGCCTCGGTG","ACAGGCTCAC","CGCCGCAGAA","GACGGTGCTC","CCGTGAAGGT","TACCTCCCTT","TAATCATCGA","CAATCTCAGT","GAGGCCCCTA","CTTTGGAATA","TGGAACTTGG","ATCTGGGGCG","AGTGCGTGAA","AGACACTAAT","CTACGGGTAT","TCAGGCCCGC","CTAGAGTCCA","TTAGATTTCC","CGCCCACCAA","ACGTTAGTCG","TAAAGAAAGT","GTGTAGCCCC","ACCAGGTTAC","TAAGGGAAAC","TCAATGTGAA","ACATTTCATT","GTGTCCATAC","GAAGACTTGC","TAGGGTACGA","TTAGCGCGGA","TAAATAAGTA","CCACTGGCAC","AGGGCCCTTG","TAGCCACTAT","GGAACCTTTG","TGTCCGATCC","CTCCGGGTAC","AAAGTGTTCG","TAGTGGGAGC","CGCAATCTTG","CTCTAATGTG","CCGCTGGATT","ACCACTTTTT","CTTGCTTGTT","TGCCTTCTAG","AAGCGATCAC","TTAGATATTA","AAGGACAATA","GTAGATTAAC","GTGTCTATCA","GATCCGACCG","GATGCACTAG","ATCGGTACCC","CTGAACTGAA","GATCGATTGA","GACTTGCTAA","TAGCCCCTAA","GCTGCTTATA","AAGCCCGGTG","TGCGCAAGAT","GTCGAGACGT","TTTTCCGACG","GATCGCGCGC","ATCGCCGCAT","CGAATGAGGT","CACTCGCATT","CATCATAAGG","CTGCTACTTA","CCGCTCTAAA","GCGGTGCCGT","GGATACCTTC","GATGTGTCTA","CGCATTTTCC","CCTGTATCGC","ATGGGACCAC","AGATCGGCGT","TCGCAGCCCT","CGGGGCCTCC","TCAGACCTGT","TCAAACTTGC","GTTCAGTTAT","GAACGTGGTA","CTCGGACAAA","GTTATTTCCG","TATTCGCTTA","TGTCCCGGGG","AGTTCCTACA","ATCCACGGAG","TAGTAGCACG","GGGTAACACA","ACCTTTATCA","CTAACATGAC","CGGTCCCTCT","CAAAATCTAT","AAAGATTCGC","CAAGAGTACC","TCATACTTTG","TCTCCTACGC","AATTAGACCC","GGAGTCATCA","GACATAGGAC","ACAGGGATAC","GTAGCTTGCG","AGCAGTTTGG","ACTGGTGGAA","GGAACTTCGG","CGACAGCATC","TCAATCCATG","CCTCCTACGA","AGAATTTCTC","ACGATTGTAG","ACTGGCTCGC","TGTATACGGG","CTTACAGAAA","GTGCTGTAGA","CGGATCCTTA","TAAAGTCCCC","CACTATTCTC","TATGGTTTTC","CTTAAAGTCC","CCCCTGGTAC","ACTACGGGTA","CACCTCGAGG","GGAAGCGCCG","TTTCTATGGA"],
	 Pos = lists:map(fun(E) -> 
			bwt_match(S,E) end,
			P),
	 lists:foreach(fun(Ps)-> io:format("~p ",[Ps]) end, Pos).

bwt_match(LastCol,Pattern) ->
	FirstCol = lists:sort(LastCol),
	LtFMap = last_to_first_map(LastCol,FirstCol),
	bwt_match(1,length(LastCol),FirstCol,LastCol,Pattern,LtFMap).

bwt_match(Top,Bottom,_FirstCol,_LastCol,[],_LtFMap) when Top =< Bottom ->
	(Bottom - Top) +1;

bwt_match(Top,Bottom,_FirstCol,_LastCol,_Pattern,_LtFMap) when Top > Bottom ->
	0;


bwt_match(Top,Bottom,FirstCol,LastCol,Pattern,LtFMap) when Top =< Bottom ->
	{Prefix,Symbol} = lists:split(length(Pattern)-1,Pattern),
	Substr = lists:sublist(LastCol,Top,Bottom-Top+1),
	case fl_pos(Substr,hd(Symbol)) of
		{FPos,LPos} ->
			NewTop = lists:nth(Top+FPos-1,LtFMap), % FPos should be position in LastCol
			NewBottom = lists:nth(Top+LPos-1,LtFMap),
			bwt_match(NewTop,NewBottom,FirstCol,LastCol,Prefix,LtFMap);
		not_found ->
			0
	end.


% Return {First,Last} position of symbol in String, none if not found.
fl_pos(String,Symbol) ->
	fl_pos1(util:idx(Symbol,String),String,Symbol).
fl_pos1(not_found,_,_) ->
	not_found;
fl_pos1(First,String,Symbol) ->
	fl_pos2(First,util:idx(Symbol,lists:reverse(String)),String,Symbol).
fl_pos2(_,not_found,_,_) ->
	not_found;
fl_pos2(First,Last,String,_Symbol) ->
	{First,(length(String)-Last)+1}.


last_to_first_map(LastCol,FirstCol) ->
	FCMap = fc_map(FirstCol),
	last_to_first_map(LastCol,FCMap,[]).

last_to_first_map([],_FCMap,Acc) -> 
	lists:reverse(Acc);

last_to_first_map([H|T],FCMap,Acc) ->
	{FCPos,Rem} = lists:split(1,proplists:get_value(H,FCMap)),
	last_to_first_map(T,[{H,Rem}|proplists:delete(H,FCMap)],[hd(FCPos)|Acc]).



fc_map(FirstCol) ->
	fc_map(FirstCol,1,[]).
fc_map([],_,Acc) ->
	Acc;
fc_map([H|T],Pos,Acc) ->
	PosL = proplists:get_value(H,Acc),
	NewAcc = case PosL of 
			undefined ->
				[{H,[Pos]}|Acc];
			L ->
				[{H,lists:sort([Pos|L])}|proplists:delete(H,Acc)]
			end,
	fc_map(T,Pos+1,NewAcc).


















