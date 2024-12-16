import numpy as np
import sys
import commpy as comm
import commpy.channelcoding.convcode as check
from wifitransmitter import WifiTransmitter




def WifiReceiver(txsignal, level):

    def level1(txsignal):

        """
        Decodes a transmitted signal to extract the original message and its length (Only Level 1).

        This function performs the following steps:
        1. Identifies extra zeros added during transmission for alignment (Note: None are added in Level 1)
        2. Decodes the length of the original message using repetitive decoding. 
        3. Deinterleaves and reconstructs the original message bits.
        4. Removes extra zeros appended to ensure message length is a multiple of 128.

        Args:
            txsignal (numpy.ndarray): The transmitted signal represented as a sequence of bits. 
                It includes:
                - Extra zeros added during transmission.
                - Encoded length bits.
                - Interleaved and repetitive coded message bits.

        Returns:
            tuple:
                - zeros_padded (int): Number of zeros padded in the transmission (always 0 for this implementation).
                - message (str): The decoded message.
                - message_length (int): The actual length of the original message.
        """

        zeros_padded = 0 #No zeros padded in level 1
        nfft = 64
        repetitive_coder_length = 3

        #Extra zeros added by zfill to make length bits of size 2*nfft
        extra_zeros = (2*nfft)%3

        length_bits_encoded = txsignal[extra_zeros:2*nfft]

        length_bits_decoded = []

        #Undo the repetitive coding
        for i in range(0,len(length_bits_encoded)-2,repetitive_coder_length):

            chunk = length_bits_encoded[i:i+repetitive_coder_length]
            num1s = np.sum(chunk)
            num0s = repetitive_coder_length- num1s

            if num1s > num0s:
                length_bits_decoded.append(1)
            else:
                length_bits_decoded.append(0)

        #Convert the list to an integer value
        length_bits_decoded = "".join(str(i) for i in length_bits_decoded)
        message_length = int(length_bits_decoded,2) #ACTUAL MESSAGE LENGTH (int)

        #Undo the Interleaving
        message_bits_encoded = txsignal[2*nfft:] 
        nsym = int(len(message_bits_encoded)/(2*nfft))

        message = ""

        for n in range(nsym):
            deinterleaved_chunk = np.reshape(np.transpose(np.reshape(message_bits_encoded[n*2*nfft:(n+1)*2*nfft],[4,-1])),[-1,])
            message_bits = (np.packbits(deinterleaved_chunk.astype(np.int8)))
            message += "".join([chr(num) for num in message_bits])

        #Getting Rid of extra 0s (These were added in transmitter to make the message length a multiple of 128)
        message = message[:message_length]

        return zeros_padded, message, message_length

    def level2():
        pass

    def level3():
        pass

    def level4():
        pass


    # return zero_padding, message, length



if __name__ == "__main__":

    #Takes in 3 args - Arg1 = message(str), Arg2 = level(int), Arg3 = SNR(int)
    txsignal = WifiTransmitter('Sup bitch? I am better than you', 1) 

    #Takes in 2 args - Arg1 = encodedmessage(np.array), Arg2 = level(int)
    WifiReceiver(txsignal=txsignal, level=1)

