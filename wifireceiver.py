import numpy as np
import sys
import commpy as comm
import commpy.channelcoding.convcode as check
from wifitransmitter import WifiTransmitter



class Demodulator:

    def __init__(self):

        #Some encoder parameters
        self.nfft = 64
        self.repetitive_coder = 3
        self.preamble = np.array([1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1,1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1])
        self.demodobj = comm.modulation.QAMModem(4)
        self.trellis = {'00':{0:{'output':'00', 'next_state':'00'},
                              1:{'output':'11', 'next_state':'01'}},
                        '01':{0:{'output':'10', 'next_state':'10'},
                              1:{'output':'01', 'next_state':'11'}},
                        '10':{0:{'output':'11', 'next_state':'00'},
                              1:{'output':'00', 'next_state':'01'}},
                        '11':{0:{'output':'01', 'next_state':'10'},
                              1:{'output':'10', 'next_state':'11'}},
                        }
        self.conv_encoder_rate = 0.5
        
    def hamm_dist(self, pair1, pair2):

        if len(pair1) != len(pair2):
            raise ValueError("Strings must be of same size")
        return sum([c1 != c2 for c1, c2 in zip(pair1,pair2)])
    
    def viterbi_decoder(self, input):
        


            




                



                







        
        

    def level1(self, txsignal):

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

        #Extra zeros added by zfill to make length bits of size 2*self.nfft
        extra_zeros = (2*self.nfft)%self.repetitive_coder

        length_bits_encoded = txsignal[extra_zeros:2*self.nfft]

        length_bits_decoded = []

        #Undo the repetitive coding
        for i in range(0,len(length_bits_encoded)-self.repetitive_coder+1,self.repetitive_coder):

            chunk = length_bits_encoded[i:i+self.repetitive_coder]
            num1s = np.sum(chunk)
            num0s = self.repetitive_coder- num1s

            if num1s > num0s:
                length_bits_decoded.append(1)
            else:
                length_bits_decoded.append(0)

        #Convert the list to an integer value
        length_bits_decoded = "".join(str(i) for i in length_bits_decoded)
        message_length = int(length_bits_decoded,2) #ACTUAL MESSAGE LENGTH (int)

        #Undo the Interleaving
        message_bits_encoded = txsignal[2*self.nfft:] 
        nsym = int(len(message_bits_encoded)/(2*self.nfft))

        message = ""

        for n in range(nsym):
            deinterleaved_chunk = np.reshape(np.transpose(np.reshape(message_bits_encoded[n*2*self.nfft:(n+1)*2*self.nfft],[4,-1])),[-1,])
            message_bits = (np.packbits(deinterleaved_chunk.astype(np.int8)))
            message += "".join([chr(num) for num in message_bits])

        #Getting Rid of extra 0s (These were added in transmitter to make the message length a multiple of 128)
        message = message[:message_length]

        return zeros_padded, message, message_length

    def level2(self, txsignal):

        #Reversing QAM Modulation
        demod_signal = self.demodobj.demodulate(txsignal, demod_type="hard")

        #Un-concatenate the preamble and length (Both of these are not conv encoded, concatenate the length back after viterbi decoding)
        # removed_preamble = demod_signal[:2*self.nfft]
        removed_length = demod_signal[2*self.nfft: 2*2*self.nfft]
        data = demod_signal[2*2*self.nfft:]

        #Start Viterbi Decoding




        
        pass

    def level3(self, txsignal):
        pass

    def level4(self, txsignal):

        pass


def WifiReceiver(txsignal, level):

    demodulator = Demodulator()

    if level == 1:
        zero_pad, message, length = demodulator.level1(txsignal=txsignal)
    elif level == 2:
        zero_pad, message, length = demodulator.level2(txsignal=txsignal)
    elif level == 3:
        zero_pad, message, length = demodulator.level3(txsignal=txsignal)
    elif level == 4:
        zero_pad, message, length = demodulator.level3(txsignal=txsignal)
    else:
        print("Invalid Level - Only Level 1 to 4 Allowed")



def main():
    #Takes in 3 args - Arg1 = message(str), Arg2 = level(int), Arg3 = SNR(int)
    txsignal = WifiTransmitter('Sup bitch? I am better than you', 2) 

    #Takes in 2 args - Arg1 = encodedmessage(np.array), Arg2 = level(int)
    WifiReceiver(txsignal=txsignal, level=2)


if __name__ == "__main__":
    main()



