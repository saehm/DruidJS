/*export class Heap {
    constructor(arr = null, accessor = (d) => d, comparator = "min") {
        this.root = null;
        this.accessor = accessor;

        if (comparator == "min") {
            this._comparator = (a, b) => a <= b;
        } else if (comparator == "max") {
            this._comparator = (a, b) => a >= b;
        } else {
            this._comparator = comparator;
        }

        if (arr && arr.length > 0) {
            let self = this;
            arr.forEach(d => self.push(d));
        }
    }

    push(element) {
        const value = this.accessor(element);
        const newNode = new Node(element, value);
        if (!this.root || this._comparator(value, this.root.value)) {
            newNode.next = this.root;
            this.root = newNode;
        } else {
            let pointer = this.root;
            while (pointer.next && !this._comparator(value, pointer.next.value)) {
                pointer = pointer.next;
            }
            newNode.next = pointer.next;
            pointer.next = newNode;
        }
        return this;
    }

    pop() {
        if (!this.root) {
            return null;
        }
        const root = this.root;
        this.root = this.root.next;
        return root;
    }

    get first() {
        return this.root;
    }

    * iterate() {
        let pointer = this.root;
        while (pointer) {
            yield pointer.element;
            pointer = pointer.next;
        }
    }

    toArray() {
        let res = [];
        let pointer = this.root;
        while (pointer) {
            res.push(pointer.element)
            pointer = pointer.next;
        }
        return res;
    }

    get length() {
        let len = 0;
        let pointer = this.root;
        while (pointer) {
            len += 1;
            pointer = pointer.next;
        }
        return len;
    }

    get empty() {
        return this.root === null;
    }
}

class Node {
    constructor(element, value) {
        this.element = element;
        this.value = value;
        this.next = null;
    }
}*/

export class Heap {
    constructor(arr = null, accessor = (d) => d, comparator = "min") {
        //this.root = null;
        this._accessor = accessor;
        this._container = []

        if (comparator == "min") {
            this._comparator = (a, b) => a <= b;
        } else if (comparator == "max") {
            this._comparator = (a, b) => a >= b;
        } else {
            this._comparator = comparator;
        }

        //console.log(arr)
        if (arr && arr.length > 0) {
            let self = this;
            arr.forEach(d => self.push(d));
        }

    }

    /* constructor(arr = null, accessor = (d) => d, comparator = "min") {
        this._accessor = accessor;

        if (comparator === "min") {
            this._comparator = (a, b) => a <= b;
        } else if (comparator === "max") {
            this._comparator = (a, b) => a >= b;
        } else {
            this._comparator = comparator;
        }

        console.log(arr)
        this._container = []
        if (arr && arr.length > 0) {
            let self = this;
            arr.forEach(d => self.push(d)) 
        }
    } */

    _get_left_child_index(index) {
        return (index * 2) + 1;
    }

    _get_right_child_index(index) {
        return (index * 2) + 2;
    }

    _get_parent_index(index) {
        return Math.floor((index - 1) / 2);
    }

    _has_parent(index) {
        return this._get_parent_index(index) >= 0;
    }

    _has_left_child(index) {
        return this._get_left_child_index(index) < this._container.length;
    }

    _has_right_child(index) {
        return this._get_right_child_index(index) < this._container.length;
    }

    _left_child(index) {
        return this._container[this._get_left_child_index(index)];
    }

    _right_child(index) {
        return this._container[this._get_right_child_index(index)];
    }

    _parent(index) {
        return this._container[this._get_parent_index(index)];
    }

    _swap(index_a, index_b) {
        //[this._container[index_b], this._container[index_a]] = [this._container[index_a], this._container[index_b]];
        let tmp = this._container[index_b];
        this._container[index_b] = this._container[index_a];
        this._container[index_a] = tmp;
    }

    push(element) {
        const value = this._accessor(element);
        const node = new Node(element, value);
        this._container.push(node);
        this._heapify_up();
        return this;
    }

    _heapify_up(start_index) {
        let index = start_index || this._container.length - 1;
        while (this._has_parent(index) && !this._comparator(this._parent(index).value, this._container[index].value)) {
            this._swap(index, this._get_parent_index(index));
            index = this._get_parent_index(index);
        }
    }

    pop() {
        if (this._container.length === 0) {
            return null;
        }
        if (this._container.length === 1) {
            return this._container.pop().element;
        }
        
        const item = this._container[0];

        this._container[0] = this._container.pop();
        this._heapify_down();

        return item;
    }

    _heapify_down(start_index=0) {
        let index = start_index;
        let next_index = null;

        while (this._has_left_child(index)) {
            if (this._has_right_child(index) && this._comparator(this._right_child(index).value, this._left_child(index).value)) {
                next_index = this._get_right_child_index(index);
            } else {
                next_index = this._get_left_child_index(index);
            }

            if (this._comparator(this._comparator(this._container[index].value, this._container[next_index].value))) {
                break;
            }

            this._swap(index, next_index);
            index = next_index;
        }
    }

    // peek
    get first() {
        return this._container.length > 0 ? this._container[0] : null;
    }

    * iterate() {
        for (let i = 0, n = this._container.length; i < n; ++i) {
            yield this._container[i].element;
        }
    }

    toArray() {
        const comparator = this._comparator;
        const accessor = this._accessor;
        let container = this._container//.map(d => d.element);
        return container.sort((a, b) => comparator(accessor(a.element), accessor(b.element)) ? -1 : 1)
        //return this._container.sort((a, b) => comparator(a.value, b.value))//.map(d => d.element);
    }

    get length() {
        return this._container.length;
    }

    get empty() {
        return this.length === 0;
    }
}

class Node {
    constructor(element, value) {
        this.element = element;
        this.value = value;
    }
}