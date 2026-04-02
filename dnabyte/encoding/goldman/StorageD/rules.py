'''
# 4bit - 2base
# 16!=20,922,789,888,000
'''
import random
import logging
log = logging.getLogger('mylog')

if 'ALL_RULES' not in locals().keys():
    ALL_RULES = list()

RULES_COUNT = 30000

base_rules_ji = {"0000": "AC", "0001": "AA", "0010": "AT", "0011": "AG", "0100": "TC", "0101": "TT", "0110": "GC",
                       "0111": "TG", "1000": "CA", "1001": "GG", "1010": "CT", "1011": "CG", "1100": "GA", "1101": "GT",
                       "1110": "CC", "1111": "TA"}
base_rules_ou = {"0000": "GG", "0001": "GC", "0010": "GT", "0011": "GA", "0100": "CT", "0101": "CC", "0110": "CG",
                       "0111": "CA", "1000": "TG", "1001": "TC", "1010": "TA", "1011": "TT", "1100": "AG", "1101": "AT",
                       "1110": "AC", "1111": "AA"}

# The fixed binary bits are incremented in the order of 0-15, and the bases are sorted (or the numbers 0-15 are sorted),
# and the corresponding relationship is obtained.
def getRules(iTotalNum=10):
    """
    Get the given number of rules
    @param iTotalNum: number of rules
    @return:
    """
    log.debug("getRules， total :{}".format(iTotalNum))
    lBase = [[x + y for x in ['A', 'T', 'C', 'G'] for y in ['G', 'C', 'T', 'A']],
             [x + y for x in ['T', 'A', 'G', 'C'] for y in ['A', 'C', 'G', 'T']],
             [x + y for x in ['G', 'A', 'C', 'T'] for y in ['C', 'A', 'G', 'T']],
             [x + y for x in ['C', 'T', 'G', 'A'] for y in ['T', 'G', 'A', 'C']]]

    lBit = [ x+y for x in ['00','01','10','11'] for y in ['00','01','10','11']]
    lResJi = getRandomRulues(16, iTotalNum, 1953)
    lResOu = getRandomRulues(16, iTotalNum, 2022)

    lPair = list()
    lPair.append((base_rules_ji, base_rules_ou))
    for index in range(iTotalNum):
        dict1,dict2 = dict(),dict()
        lRule1 = [lBase[index%4][lResJi[index][i]] for i in range(16)] #Odd Correspondence Rule
        lRule2 = [lBase[(index+2)%4][lResOu[index][i]] for i in range(16)] #even number rule

        for index,bit in enumerate(lBit):
            dict1[bit] = lRule1[index]
            dict2[bit] = lRule2[index]
        lPair.append((dict1,dict2))
    return lPair

def getRandomRulues(member_num=16, total_rules=20, seed = 2022):
    """
    Get rules by random number given random seed
    @param member_num:
    @param total_rules:
    @param seed:
    @return:
    """
    res_list = []
    template = list(range(member_num))
    random.seed(seed)
    for i in range(total_rules):
        res_list.append(random.sample(template, member_num))
    return res_list

# set rules
def setRules():
    """
    set rules in global variable
    @return:
    """
    global ALL_RULES
    if len(ALL_RULES)==0:
        ALL_RULES = getRules(RULES_COUNT)

# run the function
setRules()