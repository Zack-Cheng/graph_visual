﻿using System.Collections.Generic;

namespace GraphVisual.TreeLib
{
    public class Tree<T>
    {
        private List<TreeNode<T>> _Nodes;

        public List<TreeNode<T>> Nodes
        {
            get => _Nodes;
            set => _Nodes = value;
        }

        private TreeNode<T> _Root;

        public TreeNode<T> Root
        {
            get => _Root;
            set => _Root = value;
        }

        public Tree()
        {
            _Nodes = new List<TreeNode<T>>();
        }

        public TreeNode<T> FindNode(T w)
        {
            foreach (TreeNode<T> i in Nodes)
            {
                if (i.Value.Equals(w))
                {
                    return i;
                }
            }

            return CreateNode(w);
        }

        private TreeNode<T> CreateNode(T w)
        {
            TreeNode<T> s = new TreeNode<T>(w);
            _Nodes.Add(s);
            return s;
        }
    }
}