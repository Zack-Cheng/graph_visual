using GraphVisual.GraphD;
using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace GraphVisual.Algorithm
{
    public class GirvanNewman : IAlgorithm
    {
        private Dictionary<Edge, double> edgeBetweenness;
        private CommunityStructure Cs;
        private DGraph graph;
        private double _BestQ, Q;

        public double BestQ
        {
            get => _BestQ;
            set => _BestQ = value;
        }

        private void WriteLog(string log = "")
        {
            Debug.WriteLine(log);
        }

        public CommunityStructure FindCommunityStructure(DGraph pGraph)
        {
            graph = pGraph.Clone();

            CommunityStructure tempCS = GetCommunityStructure();

            int initCount = tempCS.Count;
            int countCommunity = initCount;

            _BestQ = 0;
            Q = 0;

            CalculateEdgeBetweenness(tempCS, graph);
            _BestQ = CalculateModularity(tempCS, pGraph);
            Cs = tempCS;
            while (true)
            {
                while (countCommunity <= initCount)
                {
                    DGraph community = RemoveMaxEdgeBetweenness(tempCS);

                    CalculateEdgeBetweenness(tempCS, community);

                    tempCS = GetCommunityStructure();
                    countCommunity = tempCS.Count;
                }

                initCount = countCommunity;

                Q = CalculateModularity(tempCS, pGraph);
                if (Q > _BestQ)
                {
                    _BestQ = Q;
                    Cs = tempCS;
                }

                if (graph.Edges.Count == 0)
                {
                    break;
                }
            }

            return Cs;
        }

        public static double CalculateModularity(CommunityStructure pCs, DGraph pOriginalGraph)
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
                    d += pOriginalGraph.FindNode(node.Label, false).AdjacencyNodes.Count;
                }

                l /= 2;

                modularity += ((double)l / numEdge) - Math.Pow(((double)d / (2 * numEdge)), 2);
            }
            return modularity;
        }

        private DGraph RemoveMaxEdgeBetweenness(CommunityStructure pTempCS)
        {
            double max = double.MinValue;
            Edge e = null;
            foreach (KeyValuePair<Edge, double> keypair in edgeBetweenness)
            {
                if (keypair.Value > max)
                {
                    max = keypair.Value;
                    e = keypair.Key;
                }
            }

            graph.Edges.Remove(e);

            e.NodeA.AdjacencyEdges.Remove(e);
            e.NodeB.AdjacencyEdges.Remove(e);

            e.NodeA.AdjacencyNodes.Remove(e.NodeB);
            e.NodeB.AdjacencyNodes.Remove(e.NodeA);

            WriteLog(" - Remove: (" + e.NodeA.Label + ", " + e.NodeB.Label + ")\t" + edgeBetweenness[e].ToString("0.00"));

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

            _CalculateEdgeBetweenness(community);
        }

        private void _CalculateEdgeBetweenness(DGraph subgraph)
        {
            if (subgraph == null)
            {
                return;
            }

            int n = subgraph.Nodes.Count;
            int MAX = 99999;

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
                    foreach (Node v in pred[w])
                    {
                        double c = ((double)(sigma[v]) / sigma[w]) * (1.0 + delta[w]);

                        Edge e = graph.FindEdge(v, w);
                        edgeBetweenness[e] += c;

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