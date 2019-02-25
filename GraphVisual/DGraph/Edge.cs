using GraphVisual.Drawing;
using System.Drawing;

namespace GraphVisual.GraphD
{
    public class Edge : IDrawable
    {
        private Node _NodeA;

        public Node NodeA
        {
            get => _NodeA;
            set => _NodeA = value;
        }

        private Node _NodeB;

        public Node NodeB
        {
            get => _NodeB;
            set => _NodeB = value;
        }

        public Edge(Node pNodeA, Node pNodeB)
        {
            _NodeA = pNodeA;
            _NodeB = pNodeB;
        }

        public void Draw(Graphics g)
        {
            g.DrawLine(Format.LinkLineColor, _NodeA.Location, _NodeB.Location);
        }
    }
}