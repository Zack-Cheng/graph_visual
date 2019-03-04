using GraphVisual.GraphD;
using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace GraphVisual.Algorithm
{
    internal class ViaSimilarity : IAlgorithm
    {
        private readonly Dictionary<Edge, double> edgeBetweenness;
        private readonly CommunityStructure Cs;
        private DGraph graph;
        private double _BestQ;
        private readonly double Q;

        public double BestQ
        {
            get => _BestQ;
            set => _BestQ = value;
        }

        private void WriteLog(string log = "")
        {
            Debug.WriteLine(log);
        }

        private struct Similarity
        {
            public Node node;
            public double similarity;
        }

        public CommunityStructure FindCommunityStructure(DGraph pGraph)
        {
            graph = pGraph.Clone();

            CommunityStructure tempCS = GetCommunityStructure();

            int initCount = tempCS.Count;
            int countCommunity = initCount;

            Dictionary<Node, List<Similarity>> hashMaxSimilarityNode = new Dictionary<Node, List<Similarity>>();

            foreach (Node node in graph.Nodes)
            {
                if (node.AdjacencyNodes.Count > 0)
                {
                    List<Similarity> lst = new List<Similarity>();

                    double max = double.MinValue;
                    foreach (Node neighborNode in node.AdjacencyNodes)
                    {
                        Similarity s = new Similarity
                        {
                            node = neighborNode,
                            similarity = ZLZSimilarity(node, neighborNode)
                        };

                        lst.Add(s);

                        if (s.similarity > max)
                        {
                            max = s.similarity;
                        }
                    }

                    int n = lst.Count;
                    for (int i = 0; i < n; i++)
                    {
                        if (lst[i].similarity < max)
                        {
                            lst.RemoveAt(i);
                            n--;
                            i--;
                        }
                    }

                    if (lst.Count >= 2)
                    {
                    }

                    hashMaxSimilarityNode.Add(node, lst);
                }

                DGraph com = new DGraph();
                com.Nodes.Add(node);
                tempCS.Add(com);
            }

            Random r = new Random();
            int nodeNumber = r.Next(0, graph.Nodes.Count);
            Node initNode = graph.Nodes[nodeNumber];

            return Cs;
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

        private double ZLZSimilarity(Node a, Node b)
        {
            if (!a.AdjacencyNodes.Contains(b))
            {
                return 0;
            }

            double z = 0f;
            List<Node> commonNeighbors = FindCommonNeighbors(a, b);

            foreach (Node c in commonNeighbors)
            {
                z += (1.0 / c.AdjacencyNodes.Count);
            }

            return z;
        }

        private List<Node> FindCommonNeighbors(Node a, Node b)
        {
            List<Node> commonNodes = new List<Node>();

            foreach (Node c in a.AdjacencyNodes)
            {
                if (!commonNodes.Contains(c) && b.AdjacencyNodes.Contains(c))
                {
                    commonNodes.Add(c);
                }
            }

            foreach (Node c in b.AdjacencyNodes)
            {
                if (!commonNodes.Contains(c) && a.AdjacencyNodes.Contains(c))
                {
                    commonNodes.Add(c);
                }
            }

            return commonNodes;
        }
    }
}