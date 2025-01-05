import numpy as np

# Function to evaluate clustering precition
def evalPrec(labels, labels_gt, K):
    N = len(labels)

    numCorrect = 0
    map = np.zeros((K, 2))
    map[:, 0] = range(K)
    for labelRes in range(K):
        resMask = (labels == labelRes)

        matched = np.zeros(K)  # i-th value of "matched" matches "i" in gt
        for labelGT in range(K):
            gtMask = (labels_gt == labelGT)
            matched[labelGT] = np.sum(resMask & gtMask)

        # find the best match
        val = np.max(matched)
        idx = np.argmax(matched)
        numCorrect += val
        # update map
        map[labelRes, 1] = idx

    # return clustering precision and label maping
    prec = numCorrect / N
    return prec, map
