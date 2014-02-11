%
% 2014-01-18 karlll <karl@ninjacontrol.com>
%

%
% Suffix tree & suffix array
% 


-module(suffix).
-compile([export_all]).


%% -------------------------------------------------------------------------- %%
%% Suffix Tree
%% -------------------------------------------------------------------------- %%


suffix_tree_test() ->
	Text = "GTGGGTGCTTGAGATTCGTCTGTTAGCATGGTATTTTACCCACTCAGATAACCAATGAAGTGCTGCTTCGGCGCCGCGAAGTGATTTTATCGAAGTGGAAGATGAATGGTGTGAGTCGCGTAGAGGTCTGGCTGCGGGGCAACGAATGGCCGATAGATATTGGACCATAGAGCGTGACATCCATGCAGGGAGTACCATATGAATCCCACCTCTTCTCGCTTCTATACTTTGTACGTATGGGGGTCTGCCGAAGCTAATATCCTACAGACAAATGCGGACCCGTACAACAATGGAAGTCACAGCGGTCCCCATTAACATCATGCTTTAGTTCACCCTTACGGGCGGAAGAGAACCGGAACGCCTGCCTGCGACGTACGCGCTCCCAGCTTCCAGGAAAAGCGCTAACCACTCAGTGGAACGGTCAGGTAAAAATAACAGACTCGCTCTTATTCAATGCCCGTGGTTGGCAACTTCTAGTCTCAGTTTACCACCGCTACCCGGAGCAGGCCAGGGGATTTAGCTCACCCTTACTGCAAACACTGAGCGCATTACGAGTATACGTTGCACAATAATAGCCGTTTGATTTAAGGAGGCGAAAAATTTACGAAAATTGATCGCAACTCGGAGTTATTTAATATTTCGTATTTATCTGGAGTAGTGGACCGAGGGATCTCAGGTGGGTTGAGTCACTGGATAGACTTCCTGTTTTCCGACGCCGACAACGGCGATGGCGCCATATGGTTTGGTTTGGAGCATTACTCTGTTAGAGACCTCTAAGCAAGCCTCGAAAAAGCACACCCCTCCCAGGGACCTTATGCAATCGAGTCTAAATCTTGAAAGCTCCTAATTCTGCATGTCAAATGTTCTATAGCTCACCATAGAGACTGTGGGTAGGCTGCACGCCGCTCATCCTCCAGGGGCCCACCGACTCGAATCAGCCTTTAGGTATAACTTCGCACTGGTGTTATCAAAAATATATGCCTCACGGAAGCCGCTAAGAAACGCCTCGTACGCGTTGAGCGCTAAACAAAGTTAGGCTTCGTCCGCTAATTCCTCTGTTCCCCGTTTTAGCTTCTTGATAAGTCGCCGGAGGAAATTCAGGTCAGTTGCGACCTGCACTAGCTCTCTTTGGATGCCGTACGATATCCTGGTACCACCGCTTTAGAAATGTTCTCCAGAGTTTATTCGAAACTTAGGATACTGGGTGTACGAAACTGGGGAAGAGTCTTGCCACCCGGAAGCAATGTCCAAAAGGCTTAAGTGGTG$",
	L = length(Text),
	Substrs = lists:zip(lists:seq(2,L),lists:duplicate(L-1,L)),
	
	{_,Edges,_} = lists:foldl(fun(E,Acc) -> {Start,End} = E, suffix_add(Acc,Text,{Start,End}) end, suffix_new(Text),Substrs),
	util:write_result_text(get_edge_labels(Edges,Text),"suffix_tree_edges.txt").


print_edge_labels(EdgeList,Text) ->
	lists:foreach(fun({_,_,{Start,End}})->
				io:format("~s~n",[string:substr(Text,Start,(End-Start)+1)])
				end, EdgeList).
				
get_edge_labels(EdgeList,Text) ->
	lists:map(fun({_,_,{Start,End}})->
				io_lib:format("~s~n",[string:substr(Text,Start,(End-Start)+1)])
				end, EdgeList).


suffix_new(Text) ->
	{
	[{1,r},{2,l,1}],
	[{1,2,{1,length(Text)}}],
	[{1,2,char_at(Text,1)}]
	}.


suffix_add(T,Text,{SStart,SEnd}=Suffix) ->
	suffix_add(T,Text,1,char_at(Text,SStart),SStart,Suffix).

suffix_add(T,Text,CurrentNode,CurrentChar,Start,Suffix) ->
	% traverse
	
	{SStart,SEnd} = Suffix,
	{NodeList,EdgeList,FirstCharList} = T,
	OutEdgesFirstChar = lists:filter(fun({N,_,C}) -> 
									case {N,C} of 
										{CurrentNode,CurrentChar} ->
											true;
										_ ->
											false
									end end,
									FirstCharList),

	case OutEdgesFirstChar of 
		[] ->
			% no edges, add new
			NewNode = {length(NodeList)+1,l,Start}, % new leaf
			NewNodeList = [NewNode|NodeList],
			NewEdge = {CurrentNode,element(1,NewNode),Suffix},
			NewEdgeList = [NewEdge|EdgeList],
			NewFirstCharList = [{CurrentNode,element(1,NewNode),CurrentChar}|FirstCharList],
			{NewNodeList,NewEdgeList,NewFirstCharList};
		
		[{CurrentNode,NodeB,CurrentChar}] ->
			
			OutEdges = lists:filter(fun(N) -> 
										case N of 
											{CurrentNode,NodeB,_NSuf} ->
												true;
											_ ->
												false
										end end,
										EdgeList),
			
			OutEdge = hd(OutEdges),
			{CurrentNode,NodeB,NSuf} = OutEdge,
			{NStart,NEnd} = NSuf,
			% longest common prefix
			LCP = lcp(Text,NSuf,Suffix),
			case LCP of % split edge or traverse?
				I when I == (NEnd-NStart)+1 -> % match whole edge, traverse to next node
					suffix_add(T,Text,NodeB,char_at(Text,SStart+LCP),Start,{SStart+LCP,SEnd});
				I -> % split edge and add node
					% create new nodes
					NewNodeA = {length(NodeList)+1,v},
					NewNodeB = {length(NodeList)+2,l,Start},
					% add to new nodelis
					NewNodeList = [NewNodeA,NewNodeB|NodeList],
					% split current out edge, add new edge and update edgelist
					[NewE1,NewE2,NewE3] = split_edge(OutEdge,LCP,{SStart+LCP,SEnd},element(1,NewNodeA),element(1,NewNodeB)),
					NewFirstCharList = update_first_char_list(Text,FirstCharList,{CurrentNode,NodeB,CurrentChar},[NewE1,NewE2,NewE3]),
					NewEdgeList = update_edge_list(EdgeList,OutEdge,[NewE1,NewE2,NewE3]),

					{NewNodeList,NewEdgeList,NewFirstCharList}
			end
	end.


% Longest Common Prefix in two substrings
% Substr = {StartPos,EndPos}

lcp(String, {Start1,End1} = S1, {Start2, End2} = S2) -> 
	lcp(String,S1,S2,0).

lcp(String, {Start1,End1} = S1, {Start2, End2} = S2, I) when Start1+I =< End1, Start2+I =< End2 -> 
	NextC1 = char_at(String,Start1+I),
	NextC2 = char_at(String,Start2+I),
	case NextC1 of 

		NextC2 ->
			lcp(String,S1,S2,I+1);

		_ ->
			I
	end;

lcp(String,S1,S2,I)  -> 
	I.

char_at(String,Pos) -> lists:nth(Pos,String).


split_edge(Edge,PrefixLen,NewSuffix,NewNodeA,NewNodeB) ->
	{A,B,S} = Edge,
	{I,J} = S,
	[{A,NewNodeA,{I,(I+PrefixLen)-1}},{NewNodeA,B,{I+PrefixLen,J}},{NewNodeA,NewNodeB,NewSuffix}].


update_first_char_list(Text,CL,OldEntry,NewEdges) ->
	NewEntries = lists:map(fun({NodeA,NodeB,{SS,_SE}}) -> {NodeA,NodeB,char_at(Text,SS)} end, NewEdges),
	NewEntries ++ lists:delete(OldEntry,CL). 

update_edge_list(EL,OldEdge,NewEdges) ->
	NewEdges ++ lists:delete(OldEdge,EL).

new_node({NodeList,_,_} = T,l,StartIdx) ->
	[{length(NodeList)+1,l,StartIdx}|NodeList].
new_node({NodeList,_,_} = T,v) ->
	[{length(NodeList)+1,v}|NodeList].




%% -------------------------------------------------------------------------- %%
%% Suffix Array
%% -------------------------------------------------------------------------- %%

% Create a suffix array from string String
suffix_array(String) ->
	suffix_array(String,0,[]).

suffix_array([],_,Acc) ->
	lists:keysort(2,Acc);
suffix_array([H|T]=S,Pos,Acc) ->
	suffix_array(T,Pos+1,[{Pos,S}|Acc]).


