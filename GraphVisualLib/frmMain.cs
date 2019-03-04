using GraphVisual.Algorithm;
using GraphVisual.Drawing;
using GraphVisual.GraphD;
using System;
using System.Drawing;
using System.Windows.Forms;

namespace GraphVisual
{
    public partial class frmMain : Form
    {
        private string graphFileName;

        public frmMain()
        {
            InitializeComponent();
            _MyDocument = new GraphDocument(drawingSpace);
        }

        private GraphDocument _MyDocument;

        public GraphDocument MyDocument
        {
            get => _MyDocument;
            set => _MyDocument = value;
        }

        public void LoadGraph(DGraph graph)
        {
            _MyDocument.Graph = graph;
            Random rand = new Random();
            for (int i = 0; i < _MyDocument.Graph.Nodes.Count; i++)
            {
                _MyDocument.Graph.Nodes[i].Location = new Point(rand.Next(Format.Setting.NodeSize, drawingSpace.Width - Format.Setting.NodeSize), rand.Next(Format.Setting.NodeSize, drawingSpace.Height - Format.Setting.NodeSize));
            }
            ReDraw();
        }

        public void LoadCS(CommunityStructure cs)
        {
            _MyDocument.IsOverlap = false;
            _MyDocument.CS = cs;

            int n = Format.Brushes.Length;
            int numCom = MyDocument.CS.Count;
            int maxCom = 3;

            int col = (numCom >= maxCom) ? 3 : numCom;
            int row = (numCom <= maxCom) ? 1 : Convert.ToInt32(Math.Ceiling((double)numCom / maxCom));

            int width = drawingSpace.Width / col;
            int height = drawingSpace.Height / row;

            Random rand = new Random();

            int margin = 60;
            for (int i = 0; i < numCom; i++)
            {
                int minWidth = width * (i % col) + margin;
                int maxWidth = minWidth + width - margin;

                int minHeight = height * (i / col) + margin;
                int maxHeight = minHeight + height - margin;

                foreach (Node node in MyDocument.CS[i].Nodes)
                {
                    Node _node = MyDocument.Graph.FindNode(node.Label);
                    _node.NodeBrush = Format.Brushes[i % n];
                    _node.Location = new Point(rand.Next(minWidth, maxWidth), rand.Next(minHeight, maxHeight));
                }
            }

            ReDraw();
        }

        private void ReDraw()
        {
            drawingSpace.Invalidate();
        }
    }
}