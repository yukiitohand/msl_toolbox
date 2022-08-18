classdef MSTBandpassList < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Red
        Green
        Blue
        Filter_L0
        Filter_L1
        Filter_L2
        Filter_L3
        Filter_L4
        Filter_L5
        Filter_L5R
        Filter_L5G
        Filter_L5B
        Filter_L6
        Filter_L6R
        Filter_L6G
        Filter_L6B
        LensL
        Filter_R0
        Filter_R1
        Filter_R2
        Filter_R3
        Filter_R4
        Filter_R4R
        Filter_R4G
        Filter_R4B
        Filter_R5
        Filter_R5R
        Filter_R5G
        Filter_R5B
        Filter_R6
        Filter_R6R
        Filter_R6G
        Filter_R6B
        LensR
        L0Red
        L0Green
        L0Blue
        R0Red
        R0Green
        R0Blue
    end
    
    methods
        function obj = MSTBandpassList(varargin)
            if isempty(varargin)
                obj.Red   = MSTQuantumEfficiency('Red');
                obj.Green = MSTQuantumEfficiency('Green');
                obj.Blue  = MSTQuantumEfficiency('Blue');

                obj.Filter_L0 = MSTBandpass('L0');
                obj.Filter_L1 = MSTBandpass('L1');
                obj.Filter_L2 = MSTBandpass('L2');
                obj.Filter_L3 = MSTBandpass('L3');
                obj.Filter_L4 = MSTBandpass('L4');
                obj.Filter_L5 = MSTBandpass('L5');
                obj.Filter_L6 = MSTBandpass('L6');
                obj.LensL = MSTBandpass('Lens','L');

                obj.Filter_R0 = MSTBandpass('R0');
                obj.Filter_R1 = MSTBandpass('R1');
                obj.Filter_R2 = MSTBandpass('R2');
                obj.Filter_R3 = MSTBandpass('R3');
                obj.Filter_R4 = MSTBandpass('R4');
                obj.Filter_R5 = MSTBandpass('R5');
                obj.Filter_R6 = MSTBandpass('R6');
                obj.LensR = MSTBandpass('Lens','R');
            else
            end
        end
        
        function [mstbp_up] = upsample(obj,varargin)
            mstbp_up       = MSTBandpassList('empty');
            
            mstbp_up.Red   = obj.Red.upsample('qe', varargin{:});
            mstbp_up.Green = obj.Green.upsample('qe', varargin{:});
            mstbp_up.Blue  = obj.Blue.upsample('qe', varargin{:});
            if ~isempty(mstbp_up.Red.responsitivity)
                mstbp_up.Red   = obj.Red.upsample('responsitivity', varargin{:});
            end
            if ~isempty(mstbp_up.Green.responsitivity)
                mstbp_up.Green   = obj.Green.upsample('responsitivity', varargin{:});
            end
            if ~isempty(mstbp_up.Blue.responsitivity)
                mstbp_up.Blue   = obj.Blue.upsample('responsitivity', varargin{:});
            end
            
            mstbp_up.Filter_L0  = obj.Filter_L0.upsample(varargin{:});
            mstbp_up.Filter_L1  = obj.Filter_L1.upsample(varargin{:});
            mstbp_up.Filter_L2  = obj.Filter_L2.upsample(varargin{:});
            mstbp_up.Filter_L3  = obj.Filter_L3.upsample(varargin{:});
            mstbp_up.Filter_L4  = obj.Filter_L4.upsample(varargin{:});
            mstbp_up.Filter_L5  = obj.Filter_L5.upsample(varargin{:});
            mstbp_up.Filter_L6  = obj.Filter_L6.upsample(varargin{:});
            mstbp_up.LensL  = obj.LensL.upsample(varargin{:});
            
            mstbp_up.Filter_R0  = obj.Filter_R0.upsample(varargin{:});
            mstbp_up.Filter_R1  = obj.Filter_R1.upsample(varargin{:});
            mstbp_up.Filter_R2  = obj.Filter_R2.upsample(varargin{:});
            mstbp_up.Filter_R3  = obj.Filter_R3.upsample(varargin{:});
            mstbp_up.Filter_R4  = obj.Filter_R4.upsample(varargin{:});
            mstbp_up.Filter_R5  = obj.Filter_R5.upsample(varargin{:});
            mstbp_up.Filter_R6  = obj.Filter_R6.upsample(varargin{:});
            mstbp_up.LensR  = obj.LensR.upsample(varargin{:});
            
        end
        
        function [] = qe2responsitivity(obj)
            obj.Red.qe2responsitivity();
            obj.Green.qe2responsitivity();
            obj.Blue.qe2responsitivity();
        end
        
        function [mstbp_sys] = get_SpectralResponse(obj,varargin)
            mstbp_sys       = MSTBandpassList('empty');
            mstbp_sys.L0Red   = MSTSpectralResponse(obj.Red,[obj.Filter_L0,obj.LensL]);
            mstbp_sys.L0Green = MSTSpectralResponse(obj.Green,[obj.Filter_L0,obj.LensL]);
            mstbp_sys.L0Blue = MSTSpectralResponse(obj.Blue,[obj.Filter_L0,obj.LensL]);
            
            mstbp_sys.Filter_L1 = MSTSpectralResponse(obj.Green,[obj.Filter_L1,obj.LensL]);
            mstbp_sys.Filter_L2 = MSTSpectralResponse(obj.Blue,[obj.Filter_L2,obj.LensL]);
            mstbp_sys.Filter_L3 = MSTSpectralResponse(obj.Red,[obj.Filter_L3,obj.LensL]);
            mstbp_sys.Filter_L4 = MSTSpectralResponse(obj.Red,[obj.Filter_L4,obj.LensL]);
            mstbp_sys.Filter_L5 = MSTSpectralResponse(obj.Green,[obj.Filter_L5,obj.LensL]);
            mstbp_sys.Filter_L6 = MSTSpectralResponse(obj.Green,[obj.Filter_L6,obj.LensL]);
            
            mstbp_sys.Filter_L5R = MSTSpectralResponse(obj.Red,[obj.Filter_L5,obj.LensL]);
            mstbp_sys.Filter_L5G = MSTSpectralResponse(obj.Green,[obj.Filter_L5,obj.LensL]);
            mstbp_sys.Filter_L5B = MSTSpectralResponse(obj.Blue,[obj.Filter_L5,obj.LensL]);
            mstbp_sys.Filter_L6R = MSTSpectralResponse(obj.Red,[obj.Filter_L6,obj.LensL]);
            mstbp_sys.Filter_L6G = MSTSpectralResponse(obj.Green,[obj.Filter_L6,obj.LensL]);
            mstbp_sys.Filter_L6B = MSTSpectralResponse(obj.Blue,[obj.Filter_L6,obj.LensL]);
            
            mstbp_sys.R0Red   = MSTSpectralResponse(obj.Red,[obj.Filter_R0,obj.LensR]);
            mstbp_sys.R0Green = MSTSpectralResponse(obj.Green,[obj.Filter_R0,obj.LensR]);
            mstbp_sys.R0Blue  = MSTSpectralResponse(obj.Blue,[obj.Filter_R0,obj.LensR]);
            mstbp_sys.Filter_R1 = MSTSpectralResponse(obj.Green,[obj.Filter_R1,obj.LensR]);
            mstbp_sys.Filter_R2 = MSTSpectralResponse(obj.Blue,[obj.Filter_R2,obj.LensR]);
            mstbp_sys.Filter_R3 = MSTSpectralResponse(obj.Red,[obj.Filter_R3,obj.LensR]);
            mstbp_sys.Filter_R4 = MSTSpectralResponse(obj.Green,[obj.Filter_R4,obj.LensR]);
            mstbp_sys.Filter_R5 = MSTSpectralResponse(obj.Green,[obj.Filter_R5,obj.LensR]);
            mstbp_sys.Filter_R6 = MSTSpectralResponse(obj.Green,[obj.Filter_R6,obj.LensR]);
            
            mstbp_sys.Filter_R4R = MSTSpectralResponse(obj.Red,[obj.Filter_R4,obj.LensR]);
            mstbp_sys.Filter_R4G = MSTSpectralResponse(obj.Green,[obj.Filter_R4,obj.LensR]);
            mstbp_sys.Filter_R4B = MSTSpectralResponse(obj.Blue,[obj.Filter_R4,obj.LensR]);
            mstbp_sys.Filter_R5R = MSTSpectralResponse(obj.Red,[obj.Filter_R5,obj.LensR]);
            mstbp_sys.Filter_R5G = MSTSpectralResponse(obj.Green,[obj.Filter_R5,obj.LensR]);
            mstbp_sys.Filter_R5B = MSTSpectralResponse(obj.Blue,[obj.Filter_R5,obj.LensR]);
            mstbp_sys.Filter_R6R = MSTSpectralResponse(obj.Red,[obj.Filter_R6,obj.LensR]);
            mstbp_sys.Filter_R6G = MSTSpectralResponse(obj.Green,[obj.Filter_R6,obj.LensR]);
            mstbp_sys.Filter_R6B = MSTSpectralResponse(obj.Blue,[obj.Filter_R6,obj.LensR]);
            
            
        end
        
        function [img_conv] = convCRISM(obj,wa,img,varargin)
            img_conv = MASTCAMImages();
            img_conv.L0Red = obj.L0Red.convolveCRISM(wa,img,varargin{:});
            img_conv.L0Blue = obj.L0Blue.convolveCRISM(wa,img,varargin{:});
            img_conv.L0Green = obj.L0Green.convolveCRISM(wa,img,varargin{:});
            
            img_conv.L1 = obj.Filter_L1.convolveCRISM(wa,img,varargin{:});
            img_conv.L2 = obj.Filter_L2.convolveCRISM(wa,img,varargin{:});
            img_conv.L3 = obj.Filter_L3.convolveCRISM(wa,img,varargin{:});
            img_conv.L4 = obj.Filter_L4.convolveCRISM(wa,img,varargin{:});
            img_conv.L5 = obj.Filter_L5.convolveCRISM(wa,img,varargin{:});
            img_conv.L6 = obj.Filter_L6.convolveCRISM(wa,img,varargin{:});
            
            img_conv.R0Red = obj.R0Red.convolveCRISM(wa,img,varargin{:});
            img_conv.R0Blue = obj.R0Blue.convolveCRISM(wa,img,varargin{:});
            img_conv.R0Green = obj.R0Green.convolveCRISM(wa,img,varargin{:});
            
            img_conv.R1 = obj.Filter_R1.convolveCRISM(wa,img,varargin{:});
            img_conv.R2 = obj.Filter_R2.convolveCRISM(wa,img,varargin{:});
            img_conv.R3 = obj.Filter_R3.convolveCRISM(wa,img,varargin{:});
            img_conv.R4 = obj.Filter_R4.convolveCRISM(wa,img,varargin{:});
            img_conv.R5 = obj.Filter_R5.convolveCRISM(wa,img,varargin{:});
            img_conv.R6 = obj.Filter_R6.convolveCRISM(wa,img,varargin{:});
            
        end
        
    end
end

