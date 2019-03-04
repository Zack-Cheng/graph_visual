using GraphVisual.Drawing;
using GraphVisual.GraphD;
using GraphVisual.TreeLib;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Drawing;
using System.Linq;

namespace GraphVisual.Algorithm
{
    public class PairBetweenness
    {
        private Dictionary<Node, double> _vertexBetweenness;

        public Dictionary<Node, double> VertexBetweenness
        {
            set => _vertexBetweenness = value;
        }

        private DGraph _Graph;

        public DGraph Graph
        {
            get => _Graph;
            set => _Graph = value;
        }

        private Node _Vertex;

        public Node Vertex
        {
            get => _Vertex;
            set => _Vertex = value;
        }

        private string _Split;

        public string Split
        {
            get => _Split;
            set => _Split = value;
        }

        private double _SplitBetweenness;

        public double SplitBetweenness
        {
            get => _SplitBetweenness;
            set => _SplitBetweenness = value;
        }

        private List<Node> _Clique;

        public List<Node> Clique
        {
            get => _Clique;
            set => _Clique = value;
        }

        private List<Edge> _Edges;

        public List<Edge> Edges
        {
            get => _Edges;
            set => _Edges = value;
        }

        public PairBetweenness()
        {
            _Clique = new List<Node>();
            _Edges = new List<Edge>();
            _Values = new Dictionary<Edge, double>();
        }

        private Dictionary<Edge, double> _Values;

        public Dictionary<Edge, double> Values
        {
            get => _Values;
            set => _Values = value;
        }

        public void SetValue(Node pNodeA, Node pNodeB, double value)
        {
            foreach (Edge e in _Edges)
            {
                if ((e.NodeA == pNodeA && e.NodeB == pNodeB) ||
                    (e.NodeA == pNodeB && e.NodeB == pNodeA))
                {
                    _Values[e] += value;
                    return;
                }
            }
        }

        public void DiviveValue(Edge e, double value)
        {
            _Values[e] /= value;
        }

        public void SetToValue(Edge e, double value)
        {
            _Values[e] = value;
        }

        private void SortPairBetweenness()
        {
            List<double> list = _Values.Values.ToList();
            list.Sort();
        }

        public void CalculateSplitBetweenness()
        {
            while (_Clique.Count > 2)
            {
                Edge e = FindEdgeMinScore();
                Node uw = new Node
                {
                    Label = e.NodeA.Label + "_" + e.NodeB.Label
                };

                _Clique.Remove(e.NodeA);
                _Clique.Remove(e.NodeB);

                _Edges.Remove(e);
                _Values.Remove(e);

                foreach (Node x in _Clique)
                {
                    Edge ux = FindEdge(e.NodeA, x);
                    Edge wx = FindEdge(e.NodeB, x);
                    Edge uwx = new Edge(uw, x);

                    _Edges.Remove(ux);
                    _Edges.Remove(wx);

                    _Edges.Add(uwx);
                    double b = _Values[ux] + _Values[wx];
                    _Values.Add(uwx, b);

                    _Values.Remove(ux);
                    _Values.Remove(wx);
                }

                _Clique.Add(uw);
            }

            if (_Edges.Count == 1)
            {
                SplitBetweenness = _Values[_Edges[0]];
            }
        }

        public Edge FindEdge(Node pNodeA, Node pNodeB)
        {
            foreach (Edge e in _Edges)
            {
                if ((e.NodeA == pNodeA && e.NodeB == pNodeB) ||
                    (e.NodeA == pNodeB && e.NodeB == pNodeA))
                {
                    return e;
                }
            }

            return null;
        }

        public Edge FindEdge(string pLabelA, string pLabelB)
        {
            foreach (Edge e in _Edges)
            {
                if ((e.NodeA.Label == pLabelA && e.NodeB.Label == pLabelB) ||
                    (e.NodeA.Label == pLabelB && e.NodeB.Label == pLabelA))
                {
                    return e;
                }
            }

            return null;
        }

        private Edge FindEdgeMinScore()
        {
            Edge _e = null;
            double min = _Values.Values.ToList().Min();
            foreach (KeyValuePair<Edge, double> keypair in _Values)
            {
                if (keypair.Value == min)
                {
                    return keypair.Key;
                    Edge e = keypair.Key;
                    _e = e;

                    string[] nodeALabel = e.NodeA.Label.Split('_');
                    string[] nodeBLabel = e.NodeB.Label.Split('_');

                    Node a = Graph.FindNode(nodeALabel[0]);
                    Node b = Graph.FindNode(nodeBLabel[0]);

                    if ((a.AdjacencyNodes.Count + b.AdjacencyNodes.Count <= Vertex.AdjacencyNodes.Count))
                    {
                        return e;
                    }
                }
            }
            return _e;
        }

        private double TestEdge(Edge pE)
        {
            Node u = pE.NodeA;
            Node w = pE.NodeB;

            double d = 0;

            foreach (Node x in _Clique)
            {
                if (x == u || x == w)
                {
                    continue;
                }

                Edge ux = FindEdge(u, x);
                Edge wx = FindEdge(w, x);

                d += _Values[ux] + _Values[wx];
            }

            return d;
        }
    }

    public class Conga : IAlgorithm
    {
        private Dictionary<Edge, double> edgeBetweenness;
        private Dictionary<Node, double> vertexBetweenness;
        private CommunityStructure Cs;
        private DGraph graph;
        private DGraph originalGraph;
        private double _BestVAD, _VAD;
        private double _BestM, _M;

        public double BestVAD => _BestVAD;

        public double BestM => _BestM;

        public DGraph Graph => graph;

        public int countSplit = 0;

        private bool _UseVAD;

        public bool UseVAD
        {
            get => _UseVAD;
            set => _UseVAD = value;
        }

        private void WriteLog(string log = "")
        {
            Debug.WriteLine(log);
        }

        public CommunityStructure FindCommunityStructure(DGraph pGraph)
        {
            originalGraph = pGraph;
            graph = pGraph.Clone();

            CommunityStructure tempCS = GetCommunityStructure();

            int initCount = tempCS.Count;
            int countCommunity = initCount;

            _BestVAD = 0;
            _VAD = 0;
            _BestM = 0;
            _M = 0;

            CalculateEdgeBetweenness(tempCS);
            CalculateVertexBetweenness(tempCS);

            WriteLog("Start CONGA algorithm");
            while (true)
            {
                while (countCommunity <= initCount)
                {
                    RemoveEdgeOrSplitVertex(tempCS);

                    tempCS = GetCommunityStructure();

                    countCommunity = tempCS.Count;

                    CalculateEdgeBetweenness(tempCS);
                    CalculateVertexBetweenness(tempCS);
                }

                initCount = countCommunity;

                if (_UseVAD)
                {
                    _VAD = CalculateVAD(tempCS);

                    if (_VAD > _BestVAD)
                    {
                        _BestVAD = _VAD;
                        Cs = tempCS;
                    }
                }
                else
                {
                    _M = CalculateModularity(tempCS, originalGraph);

                    if (_M > _BestM)
                    {
                        _BestM = _M;
                        Cs = tempCS;
                    }
                }

                if (graph.Edges.Count == 0)
                {
                    break;
                }
            }

            return Cs;
        }

        private double CalculateVAD(CommunityStructure tempCS)
        {
            double EC = 0;
            double C = 0;

            double VAD = 0;

            foreach (DGraph cs in tempCS)
            {
                double temp_EC = 0;
                foreach (Node node in cs.Nodes)
                {
                    temp_EC += node.AdjacencyEdges.Count;
                }

                EC += temp_EC;
                C += cs.Nodes.Count;
            }

            VAD = EC / C;

            return VAD;
        }

        private void RemoveEdgeOrSplitVertex(CommunityStructure tempCS)
        {
            double maxEdge = double.MinValue;
            Edge e = null;
            foreach (KeyValuePair<Edge, double> keypair in edgeBetweenness)
            {
                if (keypair.Value >= maxEdge)
                {
                    maxEdge = keypair.Value;
                    e = keypair.Key;
                }
            }

            List<Node> CandidateNode = new List<Node>();
            foreach (KeyValuePair<Node, double> keypair in vertexBetweenness)
            {
                if (keypair.Value > maxEdge && keypair.Key.AdjacencyNodes.Count >= 4)
                {
                    CandidateNode.Add(keypair.Key);
                }
            }

            if (CandidateNode.Count == 0)
            {
                RemoveEdge(e, tempCS);
            }
            else
            {
                PairBetweenness maxSplitbetweenness = CalculatePairBetweenness(CandidateNode);

                if (maxSplitbetweenness.SplitBetweenness >= maxEdge)
                {
                    SplitVertex(maxSplitbetweenness, tempCS);
                }
                else
                {
                    RemoveEdge(e, tempCS);
                }
            }
        }

        private PairBetweenness CalculatePairBetweenness(List<Node> CandidateNode)
        {
            Dictionary<Node, PairBetweenness> pair = new Dictionary<Node, PairBetweenness>();

            foreach (Node node in CandidateNode)
            {
                PairBetweenness p = new PairBetweenness
                {
                    VertexBetweenness = vertexBetweenness,
                    Graph = graph,
                    Vertex = node
                };

                for (int i = 0; i < node.AdjacencyNodes.Count; i++)
                {
                    p.Clique.Add(node.AdjacencyNodes[i]);
                    for (int j = i + 1; j < node.AdjacencyNodes.Count; j++)
                    {
                        Edge e = new Edge(node.AdjacencyNodes[i], node.AdjacencyNodes[j]);
                        p.Edges.Add(e);
                        p.Values.Add(e, 0);
                    }
                }
                pair.Add(node, p);
            }

            #region tinh_pair

            for (int i = 0; i < graph.Nodes.Count; i++)
            {
                Node a = graph.Nodes[i];
                for (int j = 0; j < graph.Nodes.Count; j++)
                {
                    if (i == j)
                    {
                        continue;
                    }

                    Node b = graph.Nodes[j];
                    List<List<TreeNode<Node>>> paths = GetAllShortestPath(graph, a, b);
                    double p = 1.0 / paths.Count;

                    foreach (List<TreeNode<Node>> path in paths)
                    {
                        for (int _i = 1; _i < path.Count - 1; _i++)
                        {
                            if (pair.ContainsKey(path[_i].Value))
                            {
                                pair[path[_i].Value].SetValue(path[_i - 1].Value, path[_i + 1].Value, p);
                            }
                        }
                    }
                }
            }

            #endregion tinh_pair

            double max = double.MinValue;
            PairBetweenness maxPair = null;
            foreach (KeyValuePair<Node, PairBetweenness> keypair in pair)
            {
                keypair.Value.CalculateSplitBetweenness();
                if (maxPair == null)
                {
                    maxPair = keypair.Value;
                    max = maxPair.SplitBetweenness;
                    continue;
                }
                if (keypair.Value.SplitBetweenness > max)
                {
                    max = keypair.Value.SplitBetweenness;
                    maxPair = keypair.Value;
                }
            }

            return maxPair;
        }

        private void GetShortestFromTree(List<List<TreeNode<Node>>> pShortestList, TreeNode<Node> node, Node a, Node b, List<TreeNode<Node>> lst = null)
        {
            List<TreeNode<Node>> lstNode = lst;
            if (node.Prev == null)
            {
                lstNode = new List<TreeNode<Node>>
                {
                    node
                };
            }

            foreach (TreeNode<Node> nextnode in node.Next)
            {
                lstNode.Add(nextnode);
                GetShortestFromTree(pShortestList, nextnode, a, b, lstNode);

                if (nextnode.Value == b)
                {
                    List<TreeNode<Node>> path = new List<TreeNode<Node>>();
                    foreach (TreeNode<Node> _node in lstNode)
                    {
                        path.Add(_node);
                    }
                    pShortestList.Add(path);
                }

                lstNode.Remove(nextnode);
            }
        }

        private List<List<TreeNode<Node>>> GetAllShortestPath(DGraph subgraph, Node a, Node b)
        {
            Queue<Node> Q = new Queue<Node>();
            Queue<TreeNode<Node>> Q_Tree = new Queue<TreeNode<Node>>();

            Dictionary<Node, int> dist = new Dictionary<Node, int>();
            int MAX = 9999;
            foreach (Node d in subgraph.Nodes)
            {
                dist.Add(d, MAX);
            }
            dist[a] = 0;
            Q.Enqueue(a);

            Tree<Node> Tree = new Tree<Node>();
            TreeNode<Node> root = new TreeNode<Node>(a);
            Tree.Nodes.Add(root);
            Tree.Root = root;

            Q_Tree.Enqueue(root);

            while (Q.Count != 0)
            {
                Node v = Q.Dequeue();
                TreeNode<Node> treeNode = Q_Tree.Dequeue();

                foreach (Node w in v.AdjacencyNodes)
                {
                    TreeNode<Node> nextTreeNode = Tree.FindNode(w);
                    if (dist[w] == MAX)
                    {
                        dist[w] = dist[v] + 1;
                        Q.Enqueue(w);

                        Q_Tree.Enqueue(nextTreeNode);
                    }
                    if (dist[w] == dist[v] + 1)
                    {
                        treeNode.Next.Add(nextTreeNode);
                        nextTreeNode.Prev = treeNode;
                    }
                }
            }

            List<List<TreeNode<Node>>> pShortestList = new List<List<TreeNode<Node>>>();
            GetShortestFromTree(pShortestList, Tree.Root, a, b);

            return pShortestList;
        }

        private DGraph RemoveEdge(Edge e, CommunityStructure pTempCS)
        {
            graph.Edges.Remove(e);

            e.NodeA.AdjacencyEdges.Remove(e);
            e.NodeB.AdjacencyEdges.Remove(e);

            e.NodeA.AdjacencyNodes.Remove(e.NodeB);
            e.NodeB.AdjacencyNodes.Remove(e.NodeA);

            WriteLog("Remove: (" + e.NodeA.Label + ", " + e.NodeB.Label + ")\t" + edgeBetweenness[e].ToString("0.00"));

            edgeBetweenness.Remove(e);

            foreach (DGraph subgraph in pTempCS)
            {
                if (subgraph.Nodes.Contains(e.NodeA))
                {
                    return subgraph;
                }
            }
            return null;
        }

        private DGraph SplitVertex(PairBetweenness pairBetweenness, CommunityStructure pTempCS)
        {
            countSplit++;
            string[] labelA = pairBetweenness.Clique[0].Label.Split('_');
            string[] labelB = pairBetweenness.Clique[1].Label.Split('_');

            Node splitNode = pairBetweenness.Vertex;

            Node splitedNodeA = new Node(splitNode.Label + "-" + string.Join("-", labelA), new Point(splitNode.Location.X + Format.Setting.NodeSize, splitNode.Location.Y + Format.Setting.NodeSize));
            Node splitedNodeB = new Node(splitNode.Label + "-" + string.Join("-", labelB), new Point(splitNode.Location.X - Format.Setting.NodeSize, splitNode.Location.Y - Format.Setting.NodeSize));

            WriteLog("Split: " + splitedNodeA.Label + "/" + splitedNodeB.Label + "\t" + pairBetweenness.SplitBetweenness.ToString("0.00"));

            graph.Nodes.Remove(splitNode);
            vertexBetweenness.Remove(splitNode);

            graph.Nodes.Add(splitedNodeA);
            graph.Nodes.Add(splitedNodeB);
            vertexBetweenness.Add(splitedNodeA, 0);
            vertexBetweenness.Add(splitedNodeB, 0);

            foreach (Edge e in splitNode.AdjacencyEdges)
            {
                Edge newe = null;
                if (labelA.Contains(e.NodeA.Label))
                {
                    newe = graph.CreateLink(e.NodeA, splitedNodeA);
                }
                else if (labelA.Contains(e.NodeB.Label))
                {
                    newe = graph.CreateLink(e.NodeB, splitedNodeA);
                }
                else if (labelB.Contains(e.NodeA.Label))
                {
                    newe = graph.CreateLink(e.NodeA, splitedNodeB);
                }
                else if (labelB.Contains(e.NodeB.Label))
                {
                    newe = graph.CreateLink(e.NodeB, splitedNodeB);
                }

                graph.Edges.Remove(e);
                edgeBetweenness.Remove(e);
                edgeBetweenness.Add(newe, 0);

                if (e.NodeA != splitNode)
                {
                    e.NodeA.AdjacencyEdges.Remove(e);
                    e.NodeA.AdjacencyNodes.Remove(splitNode);
                }
                else
                {
                    e.NodeB.AdjacencyEdges.Remove(e);
                    e.NodeB.AdjacencyNodes.Remove(splitNode);
                }
            }

            foreach (DGraph subgraph in pTempCS)
            {
                if (subgraph.Nodes.Contains(splitNode))
                {
                    return subgraph;
                }
            }
            return null;
        }

        private void CalculateVertexBetweenness(CommunityStructure tempCS, DGraph community = null)
        {
            if (vertexBetweenness == null)
            {
                vertexBetweenness = new Dictionary<Node, double>();
                foreach (Node node in graph.Nodes)
                {
                    vertexBetweenness.Add(node, 0);
                }
            }

            if (community != null)
            {
                _CalculateVertexBetweenness(community);
            }
            else
            {
                foreach (DGraph subgraph in tempCS)
                {
                    _CalculateVertexBetweenness(subgraph);
                }
            }
        }

        private void _CalculateVertexBetweenness(DGraph subgraph)
        {
            int n = subgraph.Nodes.Count;

            foreach (Node node in subgraph.Nodes)
            {
                double vb = 0;

                foreach (Edge e in node.AdjacencyEdges)
                {
                    vb += edgeBetweenness[e];
                }

                vertexBetweenness[node] = (vb / 2) - (n - 1);
            }
        }

        private double CalculateModularity(CommunityStructure pCs, DGraph pOriginalGraph)
        {
            double modularity = 0;
            int numEdge = pOriginalGraph.Edges.Count;
            foreach (DGraph csItem in pCs)
            {
                int l = 0;
                int d = 0;
                foreach (Node node in csItem.Nodes)
                {
                    l += node.AdjacencyNodes.Count;

                    if (node.Label.Contains('-'))
                    {
                        string[] temp_label = node.Label.Split('-');
                        d += pOriginalGraph.FindNode(temp_label[0], false).AdjacencyNodes.Count;
                    }
                    else
                    {
                        d += pOriginalGraph.FindNode(node.Label, false).AdjacencyNodes.Count;
                    }
                }

                l /= 2;

                modularity += ((double)l / numEdge) - Math.Pow(((double)d / (2 * numEdge)), 2);
            }
            return modularity;
        }

        private void CalculateEdgeBetweenness(CommunityStructure pCS, DGraph community = null)
        {
            if (edgeBetweenness == null)
            {
                edgeBetweenness = new Dictionary<Edge, double>();
                foreach (Edge e in graph.Edges)
                {
                    edgeBetweenness.Add(e, 0);
                }
            }

            if (community != null)
            {
                _CalculateEdgeBetweenness(community);
            }
            else
            {
                foreach (DGraph subg in pCS)
                {
                    _CalculateEdgeBetweenness(subg);
                }
            }
        }

        private void _CalculateEdgeBetweenness(DGraph subgraph)
        {
            if (subgraph == null)
            {
                return;
            }

            int n = subgraph.Nodes.Count;
            int MAX = 9999;

            foreach (Node s in subgraph.Nodes)
            {
                foreach (Edge e in s.AdjacencyEdges)
                {
                    edgeBetweenness[e] = 0;
                }
            }

            foreach (Node s in subgraph.Nodes)
            {
                Queue<Node> Q = new Queue<Node>();
                Stack<Node> S = new Stack<Node>();
                Dictionary<Node, List<Node>> pred = new Dictionary<Node, List<Node>>();

                Dictionary<Node, int> dist = new Dictionary<Node, int>();
                Dictionary<Node, int> sigma = new Dictionary<Node, int>();
                Dictionary<Node, double> delta = new Dictionary<Node, double>();

                foreach (Node d in subgraph.Nodes)
                {
                    dist.Add(d, MAX);
                    sigma.Add(d, 0);
                    delta.Add(d, 0);
                    pred.Add(d, new List<Node>());
                }

                dist[s] = 0;
                sigma[s] = 1;
                Q.Enqueue(s);

                while (Q.Count != 0)
                {
                    Node v = Q.Dequeue();
                    S.Push(v);

                    foreach (Node w in v.AdjacencyNodes)
                    {
                        if (dist[w] == MAX)
                        {
                            dist[w] = dist[v] + 1;
                            Q.Enqueue(w);
                        }
                        if (dist[w] == dist[v] + 1)
                        {
                            sigma[w] = sigma[w] + sigma[v];
                            pred[w].Add(v);
                        }
                    }
                }

                while (S.Count != 0)
                {
                    Node w = S.Pop();
                    if (s == w)
                    {
                        continue;
                    }

                    foreach (Node v in pred[w])
                    {
                        double c = ((double)(sigma[v]) / sigma[w]) * (1.0 + delta[w]);

                        Edge _e = graph.FindEdge(w, v);
                        if (_e != null)
                        {
                            edgeBetweenness[_e] += c;
                        }

                        delta[v] += c;
                    }
                }
            }
        }

        private CommunityStructure GetCommunityStructure()
        {
            CommunityStructure cs = new CommunityStructure();

            int count = 0;
            int n = graph.Nodes.Count;
            Dictionary<Node, bool> visited = new Dictionary<Node, bool>();
            for (int i = 0; i < n; i++)
            {
                visited.Add(graph.Nodes[i], false);
            }

            foreach (Node i in graph.Nodes)
            {
                if (visited[i] == false)
                {
                    count++;
                    DGraph subgraph = new DGraph();

                    Queue<Node> Q = new Queue<Node>();
                    visited[i] = true;
                    Q.Enqueue(i);

                    subgraph.Nodes.Add(i);

                    while (Q.Count != 0)
                    {
                        Node v = Q.Dequeue();
                        foreach (Node j in v.AdjacencyNodes)
                        {
                            if (visited[j] == false)
                            {
                                subgraph.Nodes.Add(j);
                                visited[j] = true;
                                Q.Enqueue(j);
                            }
                        }
                    }
                    cs.Add(subgraph);
                }
            }
            return cs;
        }
    }
}