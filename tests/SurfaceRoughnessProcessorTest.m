classdef SurfaceRoughnessProcessorTest < matlab.unittest.TestCase

    properties (Constant)
        Tol = 1e-10
    end

    methods (TestClassSetup)
        function addSourceToPath(testCase)
            src = fullfile(fileparts(mfilename('fullpath')), '..');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(src));
        end
    end

    % ── Private helpers ───────────────────────────────────────────────────────

    methods (Access = private)
        function S = make1DStruct(~, X, Z)
            S.varNames = {'X', 'Z'};
            S.obj      = struct('X', X, 'Z', Z);
            S.varProps = struct('name', {'X','Z'}, 'class', {'double','double'});
            S.type     = '1D-profile';
            S.Xdir     = 1;
            S.Ydir     = [];
            S.Lx       = max(X(:)) - min(X(:));
            S.Ly       = [];
        end

        function S = make2DStruct(~, X, Y, Z)
            S.varNames = {'X', 'Y', 'Z'};
            S.obj      = struct('X', X, 'Y', Y, 'Z', Z);
            S.varProps = struct('name', {'X','Y','Z'}, 'class', {'double','double','double'});
            S.type     = '2D-surface';
            S.Xdir     = 2;
            S.Ydir     = 1;
            S.Lx       = max(X(:)) - min(X(:));
            S.Ly       = max(Y(:)) - min(Y(:));
        end

        function a = makeSurfAnswers(~)
            a.Q1 = 'Hom';   a.Q2 = 'Irreg'; a.Q3 = 'TBL';
            a.Q4 = 'Exp';   a.Q5 = 'Sandgrain';
            a.Q6 = 'Smith'; a.Q7 = '2020';
            a.Q8 = '220Grit'; a.Q9 = '10.1000/test';
        end
    end

    % ── buildSurfStructFromMatrix ─────────────────────────────────────────────

    methods (Test)
        function testBuildSurf1DVarNames(testCase)
            A = [(0:0.1:1)', sin(2*pi*(0:0.1:1))'];
            S = SurfaceRoughnessProcessor.buildSurfStructFromMatrix(A);
            testCase.verifyEqual(S.varNames, {'X', 'Z'});
        end

        function testBuildSurf1DObjHasXandZ(testCase)
            A = [(0:0.1:1)', (0:0.1:1)'.^2];
            S = SurfaceRoughnessProcessor.buildSurfStructFromMatrix(A);
            testCase.verifyTrue(isfield(S.obj, 'X'));
            testCase.verifyTrue(isfield(S.obj, 'Z'));
        end

        function testBuildSurf2DVarNames(testCase)
            [X, Y] = meshgrid(linspace(0,1,4), linspace(0,2,3));
            Z = X .* Y;
            A = [X(:), Y(:), Z(:)];
            S = SurfaceRoughnessProcessor.buildSurfStructFromMatrix(A);
            testCase.verifyEqual(S.varNames, {'X', 'Y', 'Z'});
        end

        function testBuildSurf2DPreservesGridShape(testCase)
            Nx = 4; Ny = 3;
            [X, Y] = meshgrid(linspace(0,1,Nx), linspace(0,2,Ny));
            Z = X .* Y;
            A = [X(:), Y(:), Z(:)];
            S = SurfaceRoughnessProcessor.buildSurfStructFromMatrix(A);
            testCase.verifySize(S.obj.Z, [Ny, Nx]);
        end

        function testBuildSurfWrongColumnCountErrors(testCase)
            testCase.verifyError( ...
                @() SurfaceRoughnessProcessor.buildSurfStructFromMatrix(rand(5,4)), ...
                ?MException);
        end
    end

    % ── loadSurface ───────────────────────────────────────────────────────────

    methods (Test)
        function testLoadMatFile(testCase)
            f = [tempname '.mat'];
            testCase.addTeardown(@() delete(f));
            X = (0:0.1:1)';
            Z = sin(2*pi*X);
            save(f, 'X', 'Z');
            S = SurfaceRoughnessProcessor.loadSurface(f);
            testCase.verifyEqual(S.varNames, {'X', 'Z'});
        end

        function testLoadCsvFile(testCase)
            f = [tempname '.csv'];
            testCase.addTeardown(@() delete(f));
            writematrix([(0:0.1:1)', (0:0.1:1)'.^2], f);
            S = SurfaceRoughnessProcessor.loadSurface(f);
            testCase.verifyEqual(S.varNames, {'X', 'Z'});
        end

        function testLoadUppercaseExtension(testCase)
            f = [tempname '.CSV'];
            testCase.addTeardown(@() delete(f));
            writematrix([(0:0.1:1)', (0:0.1:1)'.^2], f);
            S = SurfaceRoughnessProcessor.loadSurface(f);
            testCase.verifyEqual(S.varNames, {'X', 'Z'});
        end

        function testLoadUnknownExtensionErrors(testCase)
            testCase.verifyError( ...
                @() SurfaceRoughnessProcessor.loadSurface('data.foobar'), ...
                ?MException);
        end
    end

    % ── determineSurfaceType ─────────────────────────────────────────────────

    methods (Test)
        function testDetects1DProfile(testCase)
            S.varNames = {'X', 'Z'};
            S.obj = struct('X', (0:0.1:1)', 'Z', zeros(11,1));
            result = SurfaceRoughnessProcessor.determineSurfaceType(S);
            testCase.verifyEqual(result.type, '1D-profile');
        end

        function testDetects2DSurface(testCase)
            [X, Y] = meshgrid(0:0.1:1, 0:0.1:1);
            S.varNames = {'X', 'Y', 'Z'};
            S.obj = struct('X', X, 'Y', Y, 'Z', zeros(size(X)));
            result = SurfaceRoughnessProcessor.determineSurfaceType(S);
            testCase.verifyEqual(result.type, '2D-surface');
        end

        function testUnknownVarCountErrors(testCase)
            S.varNames = {'Z'};
            S.obj = struct('Z', (0:0.1:1)');
            testCase.verifyError( ...
                @() SurfaceRoughnessProcessor.determineSurfaceType(S), ...
                ?MException);
        end
    end

    % ── determineXandYdir ────────────────────────────────────────────────────

    methods (Test)
        function testXdirColumnVector(testCase)
            S.varNames = {'X', 'Z'};
            S.obj = struct('X', (0:0.1:1)', 'Z', zeros(11,1));
            S.type = '1D-profile';
            result = SurfaceRoughnessProcessor.determineXandYdir(S);
            testCase.verifyEqual(result.Xdir, 1);
            testCase.verifyEmpty(result.Ydir);
        end

        function testXdirRowVector(testCase)
            S.varNames = {'X', 'Z'};
            S.obj = struct('X', 0:0.1:1, 'Z', zeros(1,11));
            S.type = '1D-profile';
            result = SurfaceRoughnessProcessor.determineXandYdir(S);
            testCase.verifyEqual(result.Xdir, 2);
            testCase.verifyEmpty(result.Ydir);
        end

        function test2DXdirVariesAlongColumns(testCase)
            % meshgrid: X(i,j)=x(j) → X varies along cols → Xdir=2, Ydir=1
            [X, Y] = meshgrid(linspace(0,1,4), linspace(0,2,3));
            S.varNames = {'X', 'Y', 'Z'};
            S.obj = struct('X', X, 'Y', Y, 'Z', zeros(size(X)));
            S.type = '2D-surface';
            result = SurfaceRoughnessProcessor.determineXandYdir(S);
            testCase.verifyEqual(result.Xdir, 2);
            testCase.verifyEqual(result.Ydir, 1);
        end
    end

    % ── roughPhysicalProp ────────────────────────────────────────────────────

    methods (Test)
        function testLxFor1DProfile(testCase)
            X = (0:0.1:1)';
            S = testCase.make1DStruct(X, zeros(11,1));
            result = SurfaceRoughnessProcessor.roughPhysicalProp(S);
            testCase.verifyEqual(result.Lx, 1.0, AbsTol=testCase.Tol);
            testCase.verifyEmpty(result.Ly);
        end

        function testLxLyFor2DSurface(testCase)
            [X, Y] = meshgrid(linspace(0,2,5), linspace(0,3,4));
            S = testCase.make2DStruct(X, Y, zeros(size(X)));
            result = SurfaceRoughnessProcessor.roughPhysicalProp(S);
            testCase.verifyEqual(result.Lx, 2.0, AbsTol=testCase.Tol);
            testCase.verifyEqual(result.Ly, 3.0, AbsTol=testCase.Tol);
        end
    end

    % ── roughnessStats ───────────────────────────────────────────────────────

    methods (Test)
        function testKminKmaxKt(testCase)
            X = linspace(0, 1, 200)';
            Z = sin(2*pi*X);
            S = testCase.make1DStruct(X, Z);
            r = SurfaceRoughnessProcessor.roughnessStats(S);
            testCase.verifyEqual(r.kmin, min(Z),          AbsTol=testCase.Tol);
            testCase.verifyEqual(r.kmax, max(Z),          AbsTol=testCase.Tol);
            testCase.verifyEqual(r.kt,   max(Z) - min(Z), AbsTol=testCase.Tol);
        end

        function testKbarKaKrmsZ(testCase)
            X = linspace(0, 1, 200)';
            Z = sin(2*pi*X);
            S = testCase.make1DStruct(X, Z);
            r = SurfaceRoughnessProcessor.roughnessStats(S);
            testCase.verifyEqual(r.kbar,   mean(Z),      AbsTol=testCase.Tol);
            testCase.verifyEqual(r.ka,     mean(abs(Z)), AbsTol=testCase.Tol);
            testCase.verifyEqual(r.krms_z, rms(Z),       AbsTol=testCase.Tol);
        end

        function testKrmsUsesDeemeanedHeight(testCase)
            X = linspace(0, 1, 200)';
            Z = sin(2*pi*X) + 5;
            S = testCase.make1DStruct(X, Z);
            r = SurfaceRoughnessProcessor.roughnessStats(S);
            testCase.verifyEqual(r.krms, rms(Z - mean(Z)), AbsTol=testCase.Tol);
        end

        function testKz5IsAverageOfFiveExtremes(testCase)
            X = linspace(0, 1, 200)';
            Z = sin(2*pi*X);
            S = testCase.make1DStruct(X, Z);
            r = SurfaceRoughnessProcessor.roughnessStats(S);
            MaxZ = sortrows(Z, 'descend');
            MinZ = sortrows(Z, 'ascend');
            expected = mean(MaxZ(1:5)) - mean(MinZ(1:5));
            testCase.verifyEqual(r.kz5, expected, AbsTol=testCase.Tol);
        end

        function testSkewnessNearZeroForSymmetricSignal(testCase)
            X = linspace(0, 1, 1000)';
            Z = sin(2*pi*X);
            S = testCase.make1DStruct(X, Z);
            r = SurfaceRoughnessProcessor.roughnessStats(S);
            testCase.verifyEqual(r.Sk, 0, AbsTol=1e-3);
        end

        function testFlatnessIsPositive(testCase)
            X = linspace(0, 1, 200)';
            Z = sin(2*pi*X);
            S = testCase.make1DStruct(X, Z);
            r = SurfaceRoughnessProcessor.roughnessStats(S);
            testCase.verifyGreaterThan(r.Fl, 0);
        end
    end

    % ── EffectiveSlope ────────────────────────────────────────────────────────

    methods (Test)
        function testFlatSurfaceZeroSlope(testCase)
            X = linspace(0, 1, 100)';
            Z = 3 * ones(100, 1);
            Es = SurfaceRoughnessProcessor.EffectiveSlope(X, Z, 1, 1.0);
            testCase.verifyEqual(Es, 0, AbsTol=testCase.Tol);
        end

        function testLinearRampKnownSlope(testCase)
            m = 2.5; N = 100;
            X = linspace(0, 1, N)';
            Z = m * X;
            L = X(end) - X(1);
            % trapz with unit spacing over N-1 values of m gives m*(N-2);
            % then * dx=1/(N-1) and / L=1 → m*(N-2)/(N-1)
            expected = m * (N-2) / (N-1);
            Es = SurfaceRoughnessProcessor.EffectiveSlope(X, Z, 1, L);
            testCase.verifyEqual(Es, expected, AbsTol=testCase.Tol);
        end
    end

    % ── CorrelationLength ─────────────────────────────────────────────────────

    methods (Test)
        function testCorrelationLengthIsPositiveFinite(testCase)
            N = 500;
            X = linspace(0, 1, N)';
            Z = sin(2*pi*10*X);
            Z = Z - mean(Z);
            Rl = SurfaceRoughnessProcessor.CorrelationLength(X, Z, 1);
            testCase.verifyTrue(isfinite(Rl));
            testCase.verifyGreaterThan(Rl, 0);
        end
    end

    % ── MeanAutoCorr_FFT ──────────────────────────────────────────────────────

    methods (Test)
        function testNormalizedPeakAtZeroLagIsOne(testCase)
            N = 64;
            A = sin(linspace(0, 4*pi, N))';
            [~, C] = SurfaceRoughnessProcessor.MeanAutoCorr_FFT(A, 1);
            center = ceil(length(C) / 2);
            testCase.verifyEqual(C(center), 1.0, AbsTol=testCase.Tol);
        end

        function testOutputLengthIs2NMinus1(testCase)
            N = 64;
            A = rand(N, 1);
            [lags, C] = SurfaceRoughnessProcessor.MeanAutoCorr_FFT(A, 1);
            testCase.verifyLength(C,    2*N - 1);
            testCase.verifyLength(lags, 2*N - 1);
        end

        function testLagsAreSymmetricAroundZero(testCase)
            N = 32;
            A = rand(N, 1);
            [lags, ~] = SurfaceRoughnessProcessor.MeanAutoCorr_FFT(A, 1);
            testCase.verifyEqual(lags(1), -(N-1), AbsTol=testCase.Tol);
            testCase.verifyEqual(lags(end),  N-1,  AbsTol=testCase.Tol);
        end
    end

    % ── findSlopeCorr ────────────────────────────────────────────────────────

    methods (Test)
        function testFindsIndexOfFirstMinimum(testCase)
            % length=9, center=5, C(5:end)=[0.2,0.1,0.1,0.2,0.3], min at index 2
            C = [0.5, 0.8, 1.0, 0.6, 0.2, 0.1, 0.1, 0.2, 0.3];
            ind = SurfaceRoughnessProcessor.findSlopeCorr(C);
            testCase.verifyEqual(ind, 2);
        end

        function testDistinctFirstMinimumPosition(testCase)
            % length=10, center=5, C(5:end)=[0.5,0.3,0.2,0.1,0.1,0.2], min at index 4
            C = [0.9, 1.0, 0.8, 0.7, 0.5, 0.3, 0.2, 0.1, 0.1, 0.2];
            ind = SurfaceRoughnessProcessor.findSlopeCorr(C);
            testCase.verifyEqual(ind, 4);
        end
    end

    % ── cleanUpStruct ────────────────────────────────────────────────────────

    methods (Test)
        function testRemovesInternalFields(testCase)
            S.obj      = struct();
            S.varProps = {};
            S.varNames = {'X', 'Z'};
            S.type     = '1D-profile';
            S.Lx       = 1.0;
            S.kmin     = -0.5;
            result = SurfaceRoughnessProcessor.cleanUpStruct(S);
            testCase.verifyFalse(isfield(result, 'obj'));
            testCase.verifyFalse(isfield(result, 'varProps'));
            testCase.verifyFalse(isfield(result, 'varNames'));
        end

        function testPreservesStatisticsFields(testCase)
            S.obj      = struct();
            S.varProps = {};
            S.varNames = {'X', 'Z'};
            S.type     = '1D-profile';
            S.Lx       = 1.0;
            S.kmin     = -0.5;
            result = SurfaceRoughnessProcessor.cleanUpStruct(S);
            testCase.verifyTrue(isfield(result, 'type'));
            testCase.verifyTrue(isfield(result, 'Lx'));
            testCase.verifyTrue(isfield(result, 'kmin'));
        end
    end

    % ── checkAnswer ──────────────────────────────────────────────────────────

    methods (Test)
        function testNonEmptyStringDoesNotError(testCase)
            SurfaceRoughnessProcessor.checkAnswer('valid');
            testCase.verifyTrue(true);
        end

        function testEmptyStringErrors(testCase)
            testCase.verifyError( ...
                @() SurfaceRoughnessProcessor.checkAnswer(''), ...
                ?MException);
        end

        function testCellWithEmptyElementErrors(testCase)
            testCase.verifyError( ...
                @() SurfaceRoughnessProcessor.checkAnswer({'valid', ''}), ...
                ?MException);
        end

        function testAllNonEmptyCellDoesNotError(testCase)
            SurfaceRoughnessProcessor.checkAnswer({'first', 'second'});
            testCase.verifyTrue(true);
        end
    end

    % ── loadQuestionnaire ────────────────────────────────────────────────────

    methods (Test)
        function testParsesAnswersSkippingCommentsAndBlanks(testCase)
            f = [tempname '.txt'];
            testCase.addTeardown(@() delete(f));
            fid = fopen(f, 'w');
            fprintf(fid, '# comment\n\nAnswer1\nAnswer2\n# another\nAnswer3\n');
            fclose(fid);
            result = SurfaceRoughnessProcessor.loadQuestionnaire(f);
            testCase.verifyEqual(result.Q1, 'Answer1');
            testCase.verifyEqual(result.Q2, 'Answer2');
            testCase.verifyEqual(result.Q3, 'Answer3');
        end

        function testFieldsNamedSequentially(testCase)
            f = [tempname '.txt'];
            testCase.addTeardown(@() delete(f));
            fid = fopen(f, 'w');
            for k = 1:5
                fprintf(fid, 'Line%d\n', k);
            end
            fclose(fid);
            result = SurfaceRoughnessProcessor.loadQuestionnaire(f);
            testCase.verifyTrue(isfield(result, 'Q1'));
            testCase.verifyTrue(isfield(result, 'Q5'));
            testCase.verifyFalse(isfield(result, 'Q6'));
        end
    end

    % ── SurfAnswers2filename ─────────────────────────────────────────────────

    methods (Test)
        function testFileNameStartsWithSurfaceStatistics(testCase)
            testCase.applyFixture(matlab.unittest.fixtures.WorkingFolderFixture);
            [~, fileName] = SurfaceRoughnessProcessor.SurfAnswers2filename(testCase.makeSurfAnswers());
            testCase.verifyTrue(startsWith(fileName, 'SurfaceStatistics_'));
        end

        function testFileNameHasXlsxExtension(testCase)
            testCase.applyFixture(matlab.unittest.fixtures.WorkingFolderFixture);
            [~, fileName] = SurfaceRoughnessProcessor.SurfAnswers2filename(testCase.makeSurfAnswers());
            testCase.verifyTrue(endsWith(fileName, '.xlsx'));
        end

        function testDirNameContainsAuthorAndYear(testCase)
            testCase.applyFixture(matlab.unittest.fixtures.WorkingFolderFixture);
            answers = testCase.makeSurfAnswers();
            [dirName, ~] = SurfaceRoughnessProcessor.SurfAnswers2filename(answers);
            testCase.verifyTrue(contains(dirName, answers.Q6));
            testCase.verifyTrue(contains(dirName, answers.Q7));
        end

        function testCreatesExpectedSubfolders(testCase)
            testCase.applyFixture(matlab.unittest.fixtures.WorkingFolderFixture);
            [dirName, ~] = SurfaceRoughnessProcessor.SurfAnswers2filename(testCase.makeSurfAnswers());
            parentDir = fileparts(dirName);
            testCase.verifyTrue(isfolder(fullfile(parentDir, 'Surfaces')));
            testCase.verifyTrue(isfolder(fullfile(parentDir, 'Flow documentation')));
            testCase.verifyTrue(isfolder(fullfile(parentDir, 'Paper')));
        end
    end

end
