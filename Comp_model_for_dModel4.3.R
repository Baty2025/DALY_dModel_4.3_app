library(DiagrammeR)

grViz("
digraph G {
  graph [rankdir=TB, fontsize=12];
  node [fontname='Helvetica', fontsize=12, penwidth=1.5];
  edge [fontname='Helvetica', fontsize=10];

  p_ncc [label='Prevalence of NCC', shape=box, style=filled, color='#000000', fillcolor='#ffffff'];
  p_n_epi [label='Prevalence of NCC with epilepsy', shape=box, style=filled, color='black', fillcolor='white'];
  p_n_hd [label='Prevalence of NCC with headache', shape=box, style=filled, color='black', fillcolor='#aad0d5'];
  m_ncc [label='NCC associated death', shape=box, style=filled, color='black', fillcolor='#e49b9b'];
  yll [label='YLLs', shape=box, style=filled, color='black', fillcolor='#e49b9b'];
  ep_tr [label='Treated NCC with epilepsy', shape=box, style=filled, color='black', fillcolor='#aad0d5'];
  ep_ntr [label='Untreated NCC with epilepsy', shape=box, style=filled, color='black', fillcolor='#aad0d5'];
  yld [label='YLDs', shape=box, style=filled, color='black', fillcolor='#aad0d5'];
  daly [label='DALYs', shape=ellipse, style=filled, color='black', fillcolor='white',penwidth=3];

  p_ncc -> p_n_epi [label='Proportion of epilepsy among NCC patients'];
  p_ncc -> p_n_hd [label='Proportion of headache among NCC patients'];
  p_n_epi -> m_ncc [label='Case-fatality of epilepsy'];
  m_ncc -> yll [label='LE'];
  yll -> daly;
  p_n_epi -> ep_tr [label='Proportion of treated epilepsy'];
  ep_tr -> yld [label='DW of treated epilepsy'];
  yld -> daly;
  p_n_epi -> ep_ntr [label='Proportion of untreated epilepsy'];
  ep_ntr -> yld [label='DW of untreated epilepsy'];
  p_n_hd -> yld [label=''];
}
")




library(DiagrammeR)

grViz("
digraph G {
  graph [rankdir=LR, fontsize=12, center=true];

  # Default node and edge styles
  node [fontname='Helvetica', fontsize=12, penwidth=1.5, style=filled];
  edge [fontname='Helvetica', fontsize=10];

  # ------------------------
  # Nodes
  # ------------------------
  p_ncc   [label='Prevalence of NCC', shape=box, fillcolor='#ffffff'];
  p_n_epi [label='Prevalence of NCC with epilepsy', shape=box, fillcolor='#cce5ff'];
  p_n_hd  [label='Prevalence of NCC with headache', shape=box, fillcolor='#cce5ff'];

  m_ncc   [label='NCC associated death', shape=box, fillcolor='#f4cccc'];
  yll     [label='YLLs', shape=box, fillcolor='#f4cccc'];

  ep_tr   [label='Treated NCC with epilepsy', shape=box, fillcolor='#d9ead3'];
  ep_ntr  [label='Untreated NCC with epilepsy', shape=box, fillcolor='#d9ead3'];
  yld     [label='YLDs', shape=box, fillcolor='#d9ead3'];

  daly    [label='DALYs', shape=ellipse, fillcolor='#ffe680', color='black', penwidth=3, fontsize=14, fontname='Helvetica-Bold', peripheries=2];

  # ------------------------
  # Edges
  # ------------------------
  p_ncc -> p_n_epi [label='Proportion of epilepsy among NCC patients', color='blue'];
  p_ncc -> p_n_hd  [label='Proportion of headache among NCC patients', color='blue'];

  p_n_epi -> m_ncc [label='Case-fatality of epilepsy', color='red'];
  m_ncc -> yll     [label='LE', color='red'];
  yll -> daly      [color='red', penwidth=2];

  p_n_epi -> ep_tr [label='Proportion treated', color='darkgreen'];
  ep_tr -> yld     [label='DW treated epilepsy', color='darkgreen'];

  p_n_epi -> ep_ntr [label='Proportion untreated', color='darkgreen'];
  ep_ntr -> yld     [label='DW untreated epilepsy', color='darkgreen'];

  p_n_hd -> yld     [label='DW headache', color='darkgreen'];
  yld -> daly       [color='darkgreen', penwidth=2];

  # ------------------------
  # Clusters (grouping)
  # ------------------------
  subgraph cluster_mortality {
    label='Mortality Pathway';
    style=rounded;
    color=lightgrey;
    m_ncc; yll;
  }

  subgraph cluster_morbidity {
    label='Morbidity Pathway';
    style=rounded;
    color=lightgrey;
    ep_tr; ep_ntr; p_n_hd; yld;
  }

  # Align YLL and YLD on the same level
  {rank=same; yll; yld}
}
")
