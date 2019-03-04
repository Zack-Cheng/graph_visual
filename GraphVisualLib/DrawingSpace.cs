using System.Drawing;
using System.Windows.Forms;

namespace GraphVisual
{
    public partial class DrawingSpace : Panel
    {
        public DrawingSpace()
        {
            SetStyle(ControlStyles.OptimizedDoubleBuffer, true);
            SetStyle(ControlStyles.UserPaint, true);
            SetStyle(ControlStyles.AllPaintingInWmPaint, true);
            Dock = DockStyle.Fill;
            BackColor = Color.White;
        }
    }
}