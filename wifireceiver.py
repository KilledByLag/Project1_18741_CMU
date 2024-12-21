import numpy as np
import sys
import commpy as comm
import commpy.channelcoding.convcode as check
from wifitransmitter import WifiTransmitter
from matplotlib import pyplot as plt



class Demodulator:
    """
    A class to handle demodulation, decoding, and signal processing for a simulated Wi-Fi receiver.
    Supports multiple decoding levels from basic message reconstruction to full OFDM frame synchronization.

    """

    def __init__(self):
        """

        Initialize the Demodulator with required parameters and objects for decoding.

        """

        #Some encoder parameters
        self.nfft = 64
        self.repetitive_coder = 3
        self.preamble = np.array([1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1,1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1])
        self.demodobj = comm.modulation.QAMModem(4)

        #Trellis Structure
        self.trellis = {'00':{'branch_1':{'input': 0,'output': -1-1j, 'prev_state':'00'},
                              'branch_2':{'input': 0,'output':  1+1j, 'prev_state':'10'}},
                        '01':{'branch_1':{'input': 1,'output':  1+1j, 'prev_state':'00'},
                              'branch_2':{'input': 1,'output': -1-1j, 'prev_state':'10'}},
                        '10':{'branch_1':{'input': 0,'output':  1-1j, 'prev_state':'01'},
                              'branch_2':{'input': 0,'output': -1+1j, 'prev_state':'11'}},
                        '11':{'branch_1':{'input': 1,'output': -1+1j, 'prev_state':'01'},
                              'branch_2':{'input': 1,'output':  1-1j, 'prev_state':'11'}},
                        }
        
    def euclidean_dist(self, num1, num2):
        """

        Compute the Euclidean distance between two complex numbers.

        Args:
            num1 (complex): First complex number.
            num2 (complex): Second complex number.

        Returns:
            float: Euclidean distance between the two numbers.

        """
        return float(np.abs(num1 - num2))
    
    def viterbi_decoder(self, input):
        """

        Perform soft Viterbi decoding on the input signal using a predefined trellis structure.

        Args:
            input (numpy.ndarray): Encoded input signal (complex values).

        Returns:
            numpy.ndarray: Decoded bit sequence.

        """
        cumulative_metrics = [{'00': {'metric': 0, 'branch': None} ,
                               '01': {'metric': np.inf,'branch': None},
                               '10': {'metric': np.inf,'branch': None},
                               '11': {'metric': np.inf,'branch': None}}]
        
        for timestep in range(1, len(input)+1, 1):
            cumulative_metrics.append({})
            for state in cumulative_metrics[timestep-1]:
                pm1 = cumulative_metrics[timestep - 1][
                self.trellis[state]['branch_1']['prev_state']]['metric'] + self.euclidean_dist(self.trellis[state]['branch_1']['output'], input[timestep-1])
                pm2 = cumulative_metrics[timestep - 1][
                self.trellis[state]['branch_2']['prev_state']]['metric'] + self.euclidean_dist(self.trellis[state]['branch_2']['output'], input[timestep-1])

                if pm1 < pm2:
                    cumulative_metrics[timestep][state] = {'metric' : pm1, 'branch' : 'branch_1'}
                else:
                    cumulative_metrics[timestep][state] = {'metric' : pm2, 'branch' : 'branch_2'}

        decoded_output = []

        #Traceback Path
        currState = None
        minMetric = np.inf

        for state in cumulative_metrics[-1]:
            if cumulative_metrics[-1][state]['metric'] <= minMetric:
                currState, minMetric = state, cumulative_metrics[-1][state]['metric']

        if currState != None:
            for timestep in range(len(input),0,-1):
                branch = cumulative_metrics[timestep][currState]['branch']
                decoded_output.append(self.trellis[currState][branch]['input'])
                currState = self.trellis[currState][branch]['prev_state']
        decoded_output = np.array(decoded_output)[::-1]

        return decoded_output

    

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

        zeros_padded_level1 = 0 #No zeros padded in level 1

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

        return zeros_padded_level1, message, message_length

    def level2(self, txsignal):
        """

        Decode the transmitted signal (Level 2) - Undo Convolutional Encoding and Pass to Level 1

        """

        zeros_padded_level2 = 0

        #Isolate Length (It has not been convolutionally coded)
        removed_length = txsignal[self.nfft:self.nfft*2]
        
        #Isolate Data
        data = txsignal[2*self.nfft: ]

        #Pass data to viterbi decoder
        data = self.viterbi_decoder(data)

        #QAM Demodulation for length
        removed_length = self.demodobj.demodulate(removed_length, 'hard')

        #Run level 1
        pass_to_level1 = np.concatenate((removed_length, data))
        zeros, message, message_length = self.level1(pass_to_level1)

        return zeros+zeros_padded_level2, message, message_length

    def level3(self, txsignal):
        """
        
        Decode the transmitted signal (Level 3) - Undo OFDM and pass to Level 2

        """

        zeros_padded_level3 = 0

        #FFTing the received signal

        output = np.zeros_like(txsignal)
        for i in range(0, len(txsignal)-self.nfft+1, 64):
            chunk = txsignal[i:i+self.nfft]
            output[i:i+self.nfft] = np.fft.fft(chunk)
        
        zeros, message, message_length = self.level2(output)

        return zeros+zeros_padded_level3, message, message_length

    def level4(self, txsignal):
        """
        
        Decode the transmitted signal (Level 3) - Perform Frame Synchronization, Remove Zeros padded and pass to Level 3

        """

        #QAM Modulating and OFDM'ing the preamble to correlate with the noisy version of the signal
        qam_modulated_preamble = self.demodobj.modulate(self.preamble.astype(bool))

        nsym_preamble = int(len(qam_modulated_preamble)/self.nfft)

        ofdm_preamble = np.zeros(len(qam_modulated_preamble), dtype=np.complex128)

        for i in range(nsym_preamble):
            chunk = qam_modulated_preamble[i*self.nfft:(i+1)*self.nfft]
            ofdm_preamble[i*self.nfft:(i+1)*self.nfft] = np.fft.ifft(chunk)
        
        #Finding Max Correlation to do frame sync
        correlation = np.correlate(ofdm_preamble, txsignal , mode='full')
        corr_max = np.argmax(correlation)

        data = txsignal[-(corr_max+1):]

        zeros_padded_level4 = len(txsignal) - len(data)

        zeros, message, message_length = self.level3(data)

        return zeros_padded_level4+zeros, message, message_length


def WifiReceiver(txsignal, level):

    demodulator = Demodulator()

    if level == 1:
        zero_pad, message, length = demodulator.level1(txsignal=txsignal)
    elif level == 2:
        zero_pad, message, length = demodulator.level2(txsignal=txsignal)
    elif level == 3:
        zero_pad, message, length = demodulator.level3(txsignal=txsignal)
    elif level == 4:
        zero_pad, message, length = demodulator.level4(txsignal=txsignal)
    else:
        print("Invalid Level - Only Level 1 to 4 Allowed")



def main():
    """
    Simulate a Wi-Fi transmitter and receiver.
    """
    zero_pad, txsignal, length = WifiTransmitter("Insert Text", 4)
    
    WifiReceiver(txsignal=txsignal, level=4)


if __name__ == "__main__":
    main()



