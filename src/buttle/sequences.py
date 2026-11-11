def length(sequence:str) -> int:
    '''
    CLI wrapper for built in len function.

    Required arguments: sequence
    '''

    return len(sequence.strip())


def find_position(sequence:str, search_string:str) -> list:
    '''
    Given a sequence and a search string, this method will find and return the position of the search string in the sequence, if it exists.
    
    Required arguments: sequence, search_string
    '''

    search_length = len(search_string)
    sequence_length = len(sequence)
    hits = []
    for i in range(0, sequence_length-search_length+1):
        sequence_part = sequence[i:i+search_length]
        if sequence_part.lower() == search_string.lower():
            hits.append(i+1)
    return hits


