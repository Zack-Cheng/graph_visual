using System.Collections.Generic;

namespace GraphVisual.TreeLib
{
    public class TreeNode<T>
    {
        private T _Value;

        public T Value
        {
            get => _Value;
            set => _Value = value;
        }

        private List<TreeNode<T>> _Nexts;

        public List<TreeNode<T>> Next
        {
            get => _Nexts;
            set => _Nexts = value;
        }

        private TreeNode<T> _Prev;

        public TreeNode<T> Prev
        {
            get => _Prev;
            set => _Prev = value;
        }

        public TreeNode(T pValue)
        {
            _Value = pValue;
            _Nexts = new List<TreeNode<T>>();
            _Prev = null;
        }
    }
}