classdef TestPhantom < matlab.unittest.TestCase
    % Unit tests for the lung phantom private helpers.

    methods (Test)
        function testConfigJsonParsing(testCase)
            jsonText = '{"H":128,"lesion":{"enable":false},"outDir":"custom_out"}';
            cfg = generateLungPhantom_pro('__private__','loadConfig', jsonText);
            testCase.verifyEqual(cfg.H, 128);
            testCase.verifyFalse(cfg.lesion.enable);
            testCase.verifyEqual(cfg.outDir, "custom_out");
        end

        function testAnimationFrameConfig(testCase)
            cfg = generateLungPhantom_pro('__private__','loadConfig', struct('frameExt','png'));
            testCase.verifyTrue(cfg.makeMP4); % legacy flag remains
            testCase.verifyEqual(cfg.frameExt, ".png");
            testCase.verifyClass(cfg.frameDir, "string");
            testCase.verifyGreaterThan(cfg.frameDigits, 0);
        end

        function testBuildLabelMapLungMask(testCase)
            cfg = generateLungPhantom_pro('__private__','loadConfig', []);
            LABEL = generateLungPhantom_pro('__private__','defaultLabels');
            [labels, lungs] = generateLungPhantom_pro('__private__','buildLabelMap', cfg, LABEL);
            testCase.verifyTrue(all(labels(lungs) == LABEL.LUNG, 'all'));
            testCase.verifyEqual(nnz(lungs & labels ~= LABEL.LUNG), 0);
        end

        function testAssignPropsDeterministic(testCase)
            cfg = generateLungPhantom_pro('__private__','loadConfig', struct());
            LABEL = generateLungPhantom_pro('__private__','defaultLabels');
            [labels, lungs] = generateLungPhantom_pro('__private__','buildLabelMap', cfg, LABEL);
            rng(cfg.seed);
            [eps1, sig1, alpha1] = generateLungPhantom_pro('__private__','assignProps', cfg, labels, lungs, LABEL);
            rng(cfg.seed);
            [eps2, sig2, alpha2] = generateLungPhantom_pro('__private__','assignProps', cfg, labels, lungs, LABEL);
            testCase.verifyEqual(eps1, eps2, 'AbsTol', 0);
            testCase.verifyEqual(sig1, sig2, 'AbsTol', 0);
            testCase.verifyEqual(alpha1, alpha2, 'AbsTol', 0);
            testCase.verifyEqual(alpha1(~lungs), zeros(sum(~lungs(:)),1), 'AbsTol', 0);
        end

        function testAssignPropsLungBlending(testCase)
            cfg = generateLungPhantom_pro('__private__','loadConfig', []);
            LABEL = generateLungPhantom_pro('__private__','defaultLabels');
            [labels, lungs] = generateLungPhantom_pro('__private__','buildLabelMap', cfg, LABEL);
            rng(cfg.seed);
            [epsr, sigma, alpha] = generateLungPhantom_pro('__private__','assignProps', cfg, labels, lungs, LABEL);
            epsVal = eps;
            testCase.verifyEqual(unique(epsr(labels == LABEL.FAT)), cfg.eps.FAT);
            testCase.verifyEqual(unique(sigma(labels == LABEL.HEART)), cfg.sig.HEART);
            testCase.verifyGreaterThanOrEqual(min(alpha(lungs), [], 'all'), cfg.alpha.min * cfg.alpha.boost - epsVal);
            testCase.verifyLessThanOrEqual(max(alpha(lungs), [], 'all'), cfg.alpha.max * cfg.alpha.boost + epsVal);
        end

        function testSmooth2NoBlur(testCase)
            img = magic(5);
            out = generateLungPhantom_pro('__private__','smooth2', img, 0);
            testCase.verifyEqual(out, img);
        end
    end
end
