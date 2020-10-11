<h1 align="center">Guidelines</h1>


<h3 align="center">Lines</h3>

```python
    [start_query_pos, start_ref_pos, end_query_pos, end_ref_pos, [dots]]
```


<h3 align="center">Actions</h3>

```python
    ["Rotation", start_line, end_line, rotation_center]

    ["Insertion", start_query_pos, start_ref_pos, insertion_length]

    ["Deletion", start_query_pos, start_ref_pos, deletion_length]

    ["Duplication", start_query_pos, start_ref_pos, duplication_length, duplication_height, before_duplication_line_index]

    ["Translocation", start_query_pos, start_ref_pos, translocation_length]
```
