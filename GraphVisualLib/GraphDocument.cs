using GraphVisual.Algorithm;
using GraphVisual.Drawing;
using GraphVisual.GraphD;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Windows.Forms;

namespace GraphVisual
{
    public class GraphDocument
    {
        private DGraph _Graph;

        public DGraph Graph
        {
            get => _Graph;
            set => _Graph = value;
        }

        private Node _HoverNode;

        public Node HoverNode
        {
            get => _HoverNode;
            set => _HoverNode = value;
        }

        private DrawingSpace _DrawControl;

        public DrawingSpace DrawControl
        {
            get => _DrawControl;
            set => _DrawControl = value;
        }

        private CommunityStructure _CS;

        public CommunityStructure CS
        {
            get => _CS;
            set => _CS = value;
        }

        private bool _IsOverlap;

        public bool IsOverlap
        {
            get => _IsOverlap;
            set => _IsOverlap = value;
        }

        private Point MouseLocation;
        private bool _IsDown;
        private bool _IsAdd;
        private bool _IsCreateLink;
        private Node NodeA;

        public bool IsCreateLink
        {
            get => _IsCreateLink;
            set => _IsCreateLink = value;
        }

        public bool IsAdd
        {
            get => _IsAdd;
            set => _IsAdd = value;
        }

        private Dictionary<Node, List<DGraph>> _OverlapNodes;

        public Dictionary<Node, List<DGraph>> OverlapNode
        {
            get => _OverlapNodes;
            set => _OverlapNodes = value;
        }

        public GraphDocument(DrawingSpace pDrawSpace)
        {
            _Graph = new DGraph();

            _DrawControl = pDrawSpace;

            _IsDown = false;
            _IsAdd = false;
            _IsCreateLink = false;

            _DrawControl.Paint += DrawingSpace_Paint;
            _DrawControl.MouseMove += DrawingSpace_MouseMove;
            _DrawControl.MouseDown += DrawingSpace_MouseDown;
            _DrawControl.MouseUp += DrawingSpace_MouseUp;
            _DrawControl.Click += DrawingSpace_Click;
        }

        private void DrawingSpace_Click(object sender, EventArgs e)
        {
            if (_IsAdd)
            {
                Node _lastNode = Graph.Nodes[Graph.Nodes.Count - 1];
                int _label = int.Parse(_lastNode.Label) + 1;
                Node _newNode = Graph.CreateNode(_label.ToString());
                _newNode.Location = new Point(MouseLocation.X, MouseLocation.Y);
                _DrawControl.Invalidate();
            }
            else if (_IsCreateLink)
            {
                if (_HoverNode == null)
                {
                    return;
                }

                if (NodeA == null)
                {
                    NodeA = _HoverNode;
                }
                else
                {
                    Graph.CreateLink(NodeA, _HoverNode);
                    NodeA = null;
                    _DrawControl.Invalidate();
                }
            }
        }

        private void DrawingSpace_MouseUp(object sender, MouseEventArgs e)
        {
            _IsDown = false;
        }

        private bool IsHover(Node pNode)
        {
            return (pNode.Location.X - Format.Setting.NodeHaftSize <= MouseLocation.X &&
                    MouseLocation.X <= pNode.Location.X + Format.Setting.NodeHaftSize &&
                    pNode.Location.Y - Format.Setting.NodeHaftSize <= MouseLocation.Y &&
                    MouseLocation.Y <= pNode.Location.Y + Format.Setting.NodeHaftSize);
        }

        private void DrawingSpace_MouseMove(object sender, System.EventArgs e)
        {
            if (_Graph == null)
            {
                return;
            }

            MouseLocation = _DrawControl.PointToClient(Cursor.Position);

            if (_IsDown == true)
            {
                _HoverNode.Location = new Point(MouseLocation.X, MouseLocation.Y);
                _DrawControl.Invalidate();
                return;
            }

            if (_HoverNode != null)
            {
                if (IsHover(_HoverNode) == false)
                {
                    _HoverNode.IsHover = false;
                    _HoverNode = null;
                    _DrawControl.Cursor = Cursors.Default;
                }
            }

            foreach (Node node in _Graph.Nodes)
            {
                if (_HoverNode != null)
                {
                    if (node == _HoverNode && _HoverNode.IsHover)
                    {
                        break;
                    }
                }

                if (IsHover(node))
                {
                    if (_HoverNode != null)
                    {
                        _HoverNode.IsHover = false;
                    }

                    _HoverNode = node;
                    _HoverNode.IsHover = true;

                    _DrawControl.Cursor = Cursors.Hand;
                }
            }

            _DrawControl.Invalidate();
        }

        private void DrawingSpace_Paint(object sender, PaintEventArgs e)
        {
            if (_Graph == null)
            {
                e.Graphics.Clear(Color.White);
                return;
            }

            e.Graphics.SmoothingMode = System.Drawing.Drawing2D.SmoothingMode.AntiAlias;

            _Graph.Draw(e.Graphics);

            if (_CS != null && _IsOverlap)
            {
                DrawOverlap(_CS, e.Graphics);
            }
        }

        private void DrawOverlap(CommunityStructure _CS, Graphics graphics)
        {
            foreach (KeyValuePair<Node, List<DGraph>> keypair in _OverlapNodes)
            {
                int n = 360 / keypair.Value.Count;
                Node node = keypair.Key;

                Rectangle bound = new Rectangle(node.Location.X - Format.Setting.NodeHaftSize, node.Location.Y - Format.Setting.NodeHaftSize, Format.Setting.NodeSize, Format.Setting.NodeSize);

                int i = 0;
                foreach (DGraph g in keypair.Value)
                {
                    graphics.FillPie(g.NodeBrush, bound, (i++ * n), n);
                }

                if (Format.Setting.ShowNodeLabel)
                {
                    graphics.DrawString(node.Label, Format.Setting.NodeLabelFont, Format.NodeLabelColor, node.Location, Format.StrFormat);
                }
            }
        }

        private void DrawingSpace_MouseDown(object sender, MouseEventArgs e)
        {
            _IsDown = (_HoverNode != null);
        }

        public void SaveToFile(string pFilename)
        {
            Bitmap bitmap = new Bitmap(_DrawControl.Width, _DrawControl.Height);
            Graphics g = Graphics.FromImage(bitmap);
            g.SmoothingMode = System.Drawing.Drawing2D.SmoothingMode.AntiAlias;
            g.FillRectangle(Brushes.White, 0, 0, _DrawControl.Width, _DrawControl.Height);
            _Graph.Draw(g);
            bitmap.Save(pFilename);
        }

        public void Dispose()
        {
            _DrawControl.Paint -= DrawingSpace_Paint;
            _DrawControl.MouseMove -= DrawingSpace_MouseMove;
            _DrawControl.MouseDown -= DrawingSpace_MouseDown;
            _DrawControl.MouseUp -= DrawingSpace_MouseUp;
            _DrawControl.Click -= DrawingSpace_Click;

            _Graph = null;
        }
    }
}