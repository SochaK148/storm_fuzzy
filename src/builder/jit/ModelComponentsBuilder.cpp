#include "src/builder/jit/ModelComponentsBuilder.h"

#include "src/models/sparse/StateLabeling.h"
#include "src/models/sparse/Dtmc.h"
#include "src/models/sparse/Ctmc.h"
#include "src/models/sparse/Mdp.h"
#include "src/models/sparse/MarkovAutomaton.h"
#include "src/models/sparse/StandardRewardModel.h"

#include "src/builder/RewardModelBuilder.h"

#include "src/settings/SettingsManager.h"
#include "src/settings/modules/CoreSettings.h"

#include "src/exceptions/WrongFormatException.h"
#include "src/utility/macros.h"
#include "src/utility/constants.h"

namespace storm {
    namespace builder {
        namespace jit {
            
            template <typename IndexType, typename ValueType>
            ModelComponentsBuilder<IndexType, ValueType>::ModelComponentsBuilder(storm::jani::ModelType const& modelType) : modelType(modelType), isDeterministicModel(storm::jani::isDeterministicModel(modelType)), isDiscreteTimeModel(storm::jani::isDiscreteTimeModel(modelType)), currentRowGroup(0), currentRow(0), transitionMatrixBuilder(std::make_unique<storm::storage::SparseMatrixBuilder<ValueType>>(0, 0, 0, true, !isDeterministicModel)) {
                
                dontFixDeadlocks = storm::settings::getModule<storm::settings::modules::CoreSettings>().isDontFixDeadlocksSet();
            }
            
            template <typename IndexType, typename ValueType>
            ModelComponentsBuilder<IndexType, ValueType>::~ModelComponentsBuilder() {
                // Intentionally left empty.
            }
            
            template <typename IndexType, typename ValueType>
            void ModelComponentsBuilder<IndexType, ValueType>::addStateBehaviour(IndexType const& stateId, StateBehaviour<IndexType, ValueType>& behaviour) {
                // Compress the choices and reduce them depending on the model type.
                behaviour.reduce(modelType);
                
                STORM_LOG_ASSERT(stateId == currentRowGroup, "Expected states in different order.");
                
                if (!isDeterministicModel) {
                    transitionMatrixBuilder->newRowGroup(currentRow);
                }
                if (!behaviour.empty()) {
                    // Add state reward entries.
                    auto stateRewardIt = behaviour.getStateRewards().begin();
                    for (auto& rewardModelBuilder : rewardModelBuilders) {
                        if (rewardModelBuilder.hasStateRewards()) {
                            rewardModelBuilder.addStateReward(*stateRewardIt);
                            ++stateRewardIt;
                        }
                    }
                    
                    for (auto const& choice : behaviour.getChoices()) {
                        // Add the elements to the transition matrix.
                        for (auto const& element : choice.getDistribution()) {
                            transitionMatrixBuilder->addNextValue(currentRow, element.getIndex(), element.getValue());
                        }
                        
                        // Add state-action reward entries.
                        auto stateActionRewardIt = choice.getRewards().begin();
                        for (auto& rewardModelBuilder : rewardModelBuilders) {
                            if (rewardModelBuilder.hasStateActionRewards()) {
                                rewardModelBuilder.addStateActionReward(*stateActionRewardIt);
                                ++stateActionRewardIt;
                            }
                        }
                        
                        // Proceed to next row.
                        ++currentRow;
                    }
                } else {
                    if (behaviour.isExpanded() && dontFixDeadlocks) {
                        STORM_LOG_THROW(false, storm::exceptions::WrongFormatException, "Error while creating sparse matrix from JANI model: found deadlock state and fixing deadlocks was explicitly disabled.");
                    } else {
                        // Add the self-loop in the transition matrix.
                        transitionMatrixBuilder->addNextValue(currentRow, currentRowGroup, storm::utility::one<ValueType>());
                        ++currentRow;
                        
                        // Add the appropriate entries for the reward models.
                        auto stateRewardIt = behaviour.getStateRewards().begin();
                        for (auto& rewardModelBuilder : rewardModelBuilders) {
                            if (rewardModelBuilder.hasStateRewards()) {
                                rewardModelBuilder.addStateReward(*stateRewardIt);
                                ++stateRewardIt;
                            }
                            if (rewardModelBuilder.hasStateActionRewards()) {
                                rewardModelBuilder.addStateActionReward(storm::utility::zero<ValueType>());
                            }
                        }
                    }
                }
                ++currentRowGroup;
            }
            
            template <typename IndexType, typename ValueType>
            storm::models::sparse::Model<ValueType, storm::models::sparse::StandardRewardModel<ValueType>>* ModelComponentsBuilder<IndexType, ValueType>::build(IndexType const& stateCount) {
                storm::storage::SparseMatrix<ValueType> transitionMatrix = this->transitionMatrixBuilder->build();
                
                // Start by building the labeling object.
                storm::models::sparse::StateLabeling stateLabeling(stateCount);
                for (auto& label : labels) {
                    stateLabeling.addLabel(label.first, std::move(label.second));
                }
                
                // Then build all reward models.
                std::unordered_map<std::string, storm::models::sparse::StandardRewardModel<ValueType>> rewardModels;
                for (auto& rewardModelBuilder : rewardModelBuilders) {
                    rewardModels.emplace(rewardModelBuilder.getName(), rewardModelBuilder.build(transitionMatrix.getRowCount(), transitionMatrix.getColumnCount(), transitionMatrix.getRowGroupCount()));
                }
                
                if (modelType == storm::jani::ModelType::DTMC) {
                    return new storm::models::sparse::Dtmc<ValueType, storm::models::sparse::StandardRewardModel<ValueType>>(std::move(transitionMatrix), std::move(stateLabeling), std::move(rewardModels));
                } else if (modelType == storm::jani::ModelType::CTMC) {
                    return new storm::models::sparse::Ctmc<ValueType, storm::models::sparse::StandardRewardModel<ValueType>>(std::move(transitionMatrix), std::move(stateLabeling), std::move(rewardModels));
                } else if (modelType == storm::jani::ModelType::MDP) {
                    return new storm::models::sparse::Mdp<ValueType, storm::models::sparse::StandardRewardModel<ValueType>>(std::move(transitionMatrix), std::move(stateLabeling), std::move(rewardModels));
                } else {
                    return nullptr;
                }
            }
            
            template <typename IndexType, typename ValueType>
            void ModelComponentsBuilder<IndexType, ValueType>::registerRewardModel(RewardModelInformation const& rewardModelInformation) {
                rewardModelBuilders.emplace_back(rewardModelInformation);
            }
            
            template <typename IndexType, typename ValueType>
            void ModelComponentsBuilder<IndexType, ValueType>::registerLabel(std::string const& name, IndexType const& stateCount) {
                labels.emplace_back(name, storm::storage::BitVector(stateCount));
            }
            
            template <typename IndexType, typename ValueType>
            void ModelComponentsBuilder<IndexType, ValueType>::addLabel(IndexType const& stateId, IndexType const& labelIndex) {
                labels[labelIndex].second.set(stateId);
            }
            
            template class ModelComponentsBuilder<uint32_t, double>;
            template class ModelComponentsBuilder<uint32_t, storm::RationalNumber>;
            template class ModelComponentsBuilder<uint32_t, storm::RationalFunction>;
            
        }
    }
}
