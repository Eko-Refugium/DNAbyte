def flatten_at_layer(nested_list, layer):
    if layer == 0:
        return nested_list
    
    flattened = []
    for element in nested_list:
        if isinstance(element, list) and layer > 0:
            flattened.extend(flatten_at_layer(element, layer - 1))
        else:
            flattened.append(element)
    
    return flattened