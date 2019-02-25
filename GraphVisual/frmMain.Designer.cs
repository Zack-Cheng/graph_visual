namespace GraphVisual
{
    partial class frmMain
    {
        /// <summary>
        /// Required designer variable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows Form Designer generated code

        /// <summary>
        /// Required method for Designer support - do not modify
        /// the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            this.drawingSpace = new GraphVisual.DrawingSpace();
            this.SuspendLayout();
            // 
            // drawingSpace
            // 
            this.drawingSpace.BackColor = System.Drawing.Color.White;
            this.drawingSpace.Dock = System.Windows.Forms.DockStyle.Fill;
            this.drawingSpace.Location = new System.Drawing.Point(0, 0);
            this.drawingSpace.Name = "drawingSpace";
            this.drawingSpace.Size = new System.Drawing.Size(1146, 876);
            this.drawingSpace.TabIndex = 0;
            // 
            // frmMain
            // 
            this.ClientSize = new System.Drawing.Size(1146, 876);
            this.Controls.Add(this.drawingSpace);
            this.FormBorderStyle = System.Windows.Forms.FormBorderStyle.FixedToolWindow;
            this.Name = "frmMain";
            this.StartPosition = System.Windows.Forms.FormStartPosition.CenterParent;
            this.ResumeLayout(false);

        }

        #endregion
        private System.Windows.Forms.ToolStripMenuItem addNodeToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem deleteNodeToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem createLinkToolStripMenuItem;
        private DrawingSpace drawingSpace;
    }
}

