"""Merges overlapping intervals

both overlaps and inclusions count
"""


def merge_overlapping(ints):  # ints - list of intervals
    """Merges overlapping intervals

    Args:
        list of tuples, each tuple contains 2 integers,
        representing the start and the end of the interval
    Returns:
        list of tuples, representing merged intervals

    Both overlaps and inclusions considered

    """
    ints.sort()  # sorting enables merging in 1 attempt
    for i in range(0, len(ints), 1):
        for j in range(i+1, len(ints), 1):  # starting from 0 is not necessary
            set1 = set(range(ints[i][0], ints[i][1] + 1, 1))
            set2 = set(range(ints[j][0], ints[j][1] + 1, 1))
            if not set(set1).isdisjoint(set(set2)):  # overlap found
                start = min(set1 | set2)  # union
                end = max(set1 | set2)    # union
                ints[i] = (start, end)  # merging
                ints[j] = (start, end)  # merging; 2nd copy will be annihilated later
    # print(ints)  # for comparison with initial list - before annihilation
    return list(sorted(set(ints)))  # annihilation


def _main():  # example of usage
    are_coord_lst_str = [
        '102648798:102648804', '102648867:102648873', '102648873:102648878',
        '102648932:102648937', '102648943:102648948', '102649036:102649041',
        '102649106:102649111', '102649107:102649112', '102649116:102649125',
        '102649117:102649124', '102649118:102649123', '102649315:102649322',
        '102649394:102649399', '102649452:102649459', '102649456:102649461',
        '102649457:102649462', '102649464:102649469', '102649466:102649473',
        '102649493:102649504', '102649494:102649503', '102649495:102649502',
        '102649496:102649501', '102649573:102649578', '102649573:102649580',
        '102649574:102649579', '102649575:102649580', '102649596:102649602',
        '102649659:102649664', '102649686:102649691', '102649723:102649729',
        '102649838:102649843', '102649878:102649883', '102649878:102649885',
        '102649878:102649887', '102649879:102649884', '102649879:102649886',
        '102649880:102649885', '102649896:102649901']

    are_coord_lst_tpl = [tuple(map(int, interval.split(':'))) for interval in are_coord_lst_str]
    # pairs of integers are better to process than strings - that's why the function takes list of tuples
    print(are_coord_lst_tpl)  # for comparison(1)
    are_coord_lst_tpl = merge_overlapping(are_coord_lst_tpl)
    print(are_coord_lst_tpl)  # for comparison(2)


if __name__ == "__main__":
    _main()
