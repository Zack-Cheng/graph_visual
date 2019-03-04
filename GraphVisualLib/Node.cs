using GraphVisual.Drawing;
using System.Collections.Generic;
using System.Drawing;

namespace GraphVisual.GraphD
{
    public class Node : IDrawable
    {
        private string _Label;

        public string Label
        {
            get => _Label;
            set => _Label = value;
        }

        private Point _Location;

        public Point Location
        {
            get => _Location;
            set => _Location = value;
        }

        private List<Node> _AdjacencyNodes;

        public List<Node> AdjacencyNodes
        {
            get => _AdjacencyNodes;
            set => _AdjacencyNodes = value;
        }

        private List<Edge> _AdjacencyEdges;

        public List<Edge> AdjacencyEdges
        {
            get => _AdjacencyEdges;
            set => _AdjacencyEdges = value;
        }

        private bool _IsHover;

        public bool IsHover
        {
            get => _IsHover;
            set => _IsHover = value;
        }

        private Brush _NodeBrush;

        public Brush NodeBrush
        {
            get => _NodeBrush;
            set => _NodeBrush = value;
        }

        public Node(string pLabel, Point pLocation)
        {
            _AdjacencyNodes = new List<Node>();
            _AdjacencyEdges = new List<Edge>();
            _Label = pLabel;
            _Location = pLocation;
            _IsHover = false;
        }

        public Node()
        {
            _AdjacencyNodes = new List<Node>();
            _AdjacencyEdges = new List<Edge>();
            _IsHover = false;
            _NodeBrush = Format.NodeBackground;
        }

        public void Draw(Graphics g)
        {
            Rectangle bound = new Rectangle(_Location.X - Format.Setting.NodeHaftSize, _Location.Y - Format.Setting.NodeHaftSize, Format.Setting.NodeSize, Format.Setting.NodeSize);

            g.DrawEllipse(Pens.Black, bound);
            if (_IsHover == false)
            {
                g.FillEllipse(_NodeBrush, bound);
            }
            else
            {
                g.FillEllipse(Format.NodeHoverBackground, bound);
            }

            if (Format.Setting.ShowNodeLabel)
            {
                g.DrawString(_Label, Format.Setting.NodeLabelFont, Format.NodeLabelColor, _Location, Format.StrFormat);
            }
        }
    }
}